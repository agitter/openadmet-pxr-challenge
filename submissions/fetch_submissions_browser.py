#!/usr/bin/env python3
"""
fetch_submissions_browser.py

Download participant method-summary files from the OpenADMET PXR challenge,
handling the diverse link formats (GitHub, gists, Google Docs, HF, OneDrive,
PDFs, blogs) using a real Chromium browser via Playwright.

Strategy per URL type:
  - GitHub blob / raw / gist  -> fetch raw text, save as .md/.txt/.py
  - .pdf links                -> download the PDF bytes directly
  - Google Docs               -> export to PDF via the public export URL if
                                 possible; else render page to PDF
  - HTML pages / blogs / HF   -> render in Chromium, save page as PDF + HTML
  - Anything requiring login  -> logged to manual_todo.txt (cannot bypass)

Setup (one time):
    pip install playwright pandas
    playwright install chromium

Usage:
    python fetch_submissions_browser.py \
        --csv 2026-06-30.csv \
        --outdir submissions_downloaded

Output:
    submissions_downloaded/<username>.<ext>   one file per accessible entry
    submissions_downloaded/_manifest.csv      what was fetched / how / status
    submissions_downloaded/manual_todo.txt    URLs you must open yourself
"""
# Not run, may revisit later

import argparse
import re
import time
from pathlib import Path
from urllib.parse import urlparse, parse_qs

import pandas as pd


def classify(url):
    u = urlparse(url)
    host = u.netloc.lower()
    path = u.path.lower()
    if "github.com" in host and "/blob/" in path:
        return "github_blob"
    if "raw.githubusercontent.com" in host:
        return "github_raw"
    if "gist.github.com" in host:
        return "gist"
    if "github.com" in host:
        return "github_repo"
    if path.endswith(".pdf"):
        return "pdf"
    if "docs.google.com" in host:
        return "gdoc"
    if "drive.google.com" in host:
        return "gdrive"
    if "1drv.ms" in host or "onedrive" in host:
        return "onedrive"
    if "huggingface.co" in host or "hf.co" in host:
        return "hf"
    return "html"


def github_blob_to_raw(url):
    # https://github.com/u/r/blob/main/f.md -> raw.githubusercontent.com/u/r/main/f.md
    return (url.replace("https://github.com/", "https://raw.githubusercontent.com/")
               .replace("/blob/", "/"))


def gdoc_export_pdf(url):
    """Turn a Google Docs edit/view URL into its PDF export URL, if it has
    a document id. Works only for publicly-shared docs."""
    m = re.search(r"/document/d/([a-zA-Z0-9_-]+)", url)
    if m:
        return f"https://docs.google.com/document/d/{m.group(1)}/export?format=pdf"
    m = re.search(r"/presentation/d/([a-zA-Z0-9_-]+)", url)
    if m:
        return f"https://docs.google.com/presentation/d/{m.group(1)}/export/pdf"
    return None


def safe_name(username, suffix):
    base = re.sub(r"[^\w.-]", "_", str(username))
    return f"{base}{suffix}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", required=True)
    ap.add_argument("--outdir", default="submissions_downloaded")
    ap.add_argument("--timeout", type=int, default=30000,
                    help="Per-page timeout in ms")
    ap.add_argument("--headless", action="store_true", default=True)
    args = ap.parse_args()

    from playwright.sync_api import sync_playwright

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.csv)
    sub = df[df["status"] == "submitted"].dropna(subset=["url"]).copy()
    print(f"{len(sub)} submitted entries with URLs")

    manifest = []
    manual = []

    with sync_playwright() as p:
        browser = p.chromium.launch(headless=args.headless)
        ctx = browser.new_context(accept_downloads=True,
                                   user_agent=("Mozilla/5.0 (X11; Linux x86_64) "
                                               "AppleWebKit/537.36 (KHTML, like "
                                               "Gecko) Chrome/120 Safari/537.36"))
        page = ctx.new_page()
        page.set_default_timeout(args.timeout)

        for _, r in sub.iterrows():
            user, url = r["username"], r["url"]
            kind = classify(url)
            status, saved = "pending", None
            try:
                if kind in ("github_blob", "github_raw", "gist"):
                    raw = (github_blob_to_raw(url) if kind == "github_blob"
                           else url)
                    if kind == "gist":
                        raw = url.rstrip("/") + "/raw"
                    resp = ctx.request.get(raw)
                    if resp.ok:
                        text = resp.text()
                        ext = Path(urlparse(raw).path).suffix or ".txt"
                        fp = outdir / safe_name(user, ext)
                        fp.write_text(text, errors="replace")
                        saved, status = fp.name, "ok"
                    else:
                        status = f"http_{resp.status}"

                elif kind == "pdf":
                    resp = ctx.request.get(url)
                    if resp.ok:
                        fp = outdir / safe_name(user, ".pdf")
                        fp.write_bytes(resp.body())
                        saved, status = fp.name, "ok"
                    else:
                        status = f"http_{resp.status}"

                elif kind == "github_repo":
                    # Try common report filenames on the default branch
                    base = url.rstrip("/")
                    m = re.search(r"github\.com/([^/]+)/([^/]+)", base)
                    got = False
                    if m:
                        u_, r_ = m.group(1), m.group(2)
                        for branch in ("main", "master"):
                            for fname in ("README.md", "METHOD_REPORT.md",
                                          "WRITEUP.md", "report.md",
                                          "METHODS.md", "method.md"):
                                raw = (f"https://raw.githubusercontent.com/"
                                       f"{u_}/{r_}/{branch}/{fname}")
                                rr = ctx.request.get(raw)
                                if rr.ok:
                                    fp = outdir / safe_name(
                                        user, f"_{fname}")
                                    fp.write_text(rr.text(), errors="replace")
                                    saved, status, got = fp.name, "ok", True
                                    break
                            if got:
                                break
                    if not got:
                        # fall back to rendering the repo page to PDF
                        page.goto(base, wait_until="networkidle")
                        fp = outdir / safe_name(user, "_repo.pdf")
                        page.pdf(path=str(fp))
                        saved, status = fp.name, "ok_rendered"

                elif kind in ("gdoc",):
                    export = gdoc_export_pdf(url)
                    done = False
                    if export:
                        resp = ctx.request.get(export)
                        # public docs return a PDF; private return an HTML login
                        ctype = resp.headers.get("content-type", "")
                        if resp.ok and "pdf" in ctype:
                            fp = outdir / safe_name(user, ".pdf")
                            fp.write_bytes(resp.body())
                            saved, status, done = fp.name, "ok_export", True
                    if not done:
                        status = "needs_login (private gdoc)"
                        manual.append((user, url, status))

                elif kind in ("gdrive", "onedrive"):
                    status = "needs_login (drive/onedrive)"
                    manual.append((user, url, status))

                else:  # html, hf, blogs
                    page.goto(url, wait_until="networkidle")
                    time.sleep(1.5)  # let late JS settle
                    fp_pdf = outdir / safe_name(user, ".pdf")
                    page.pdf(path=str(fp_pdf))
                    fp_html = outdir / safe_name(user, ".html")
                    fp_html.write_text(page.content(), errors="replace")
                    saved, status = fp_pdf.name, "ok_rendered"

            except Exception as e:
                status = f"error: {str(e)[:80]}"

            manifest.append({"username": user, "url": url, "kind": kind,
                             "saved_as": saved, "status": status})
            print(f"  {user:22s} {kind:12s} {status}")

        browser.close()

    pd.DataFrame(manifest).to_csv(outdir / "_manifest.csv", index=False)
    if manual:
        with open(outdir / "manual_todo.txt", "w") as f:
            for user, url, why in manual:
                f.write(f"{user}\t{why}\t{url}\n")
    ok = sum(1 for m in manifest if m["status"].startswith("ok"))
    print(f"\nDone. {ok}/{len(manifest)} fetched. "
          f"{len(manual)} need manual retrieval (see manual_todo.txt).")
    print(f"Files + _manifest.csv in {outdir}/")


if __name__ == "__main__":
    main()
