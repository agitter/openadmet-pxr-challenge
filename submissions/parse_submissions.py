#!/usr/bin/env python3
"""
parse_submissions.py

Extract participant method-summary links from the OpenADMET PXR challenge
Gradio /config JSON (saved to a local file).

The config is large; we don't parse it as strict JSON. Instead we scan for
the [username, timestamp, link_html] triples that appear in the leaderboard
dataframe component, and pull out the href URL from each link cell.

Usage:
    # Save the /config page source to a file first, e.g. config.json, then:
    python parse_submissions.py --config config.json --out submissions.csv
"""

import argparse
import csv
import json
import re
from pathlib import Path


HREF_RE = re.compile(r'href=[\\"\']*([^\\"\'>\s]+)')


def extract_from_json(path):
    """Try strict JSON first: walk the structure looking for lists of
    3-element rows whose 3rd element is a link cell or 'Not submitted'."""
    text = Path(path).read_text(errors="replace")
    rows = []

    # Attempt structured parse; fall back to regex if it fails.
    try:
        data = json.loads(text)
        rows = walk_for_rows(data)
        if rows:
            return rows
    except Exception:
        pass

    # Regex fallback: find [ "name", "timestamp UTC", "cell" ] triples.
    # Timestamps look like 2026-06-30 05:37 UTC.
    triple_re = re.compile(
        r'\[\s*"([^"]+)"\s*,\s*"([^"]*UTC)"\s*,\s*"((?:[^"\\]|\\.)*)"\s*\]')
    for m in triple_re.finditer(text):
        user, ts, cell = m.group(1), m.group(2), m.group(3)
        rows.append(parse_row(user, ts, cell))
    return rows


def walk_for_rows(obj):
    """Recursively search a parsed JSON structure for leaderboard rows."""
    found = []
    if isinstance(obj, list):
        # Is this a row: [str, str-with-UTC, str]?
        if (len(obj) == 3 and all(isinstance(x, str) for x in obj)
                and "UTC" in obj[1]):
            found.append(parse_row(obj[0], obj[1], obj[2]))
        else:
            for item in obj:
                found.extend(walk_for_rows(item))
    elif isinstance(obj, dict):
        for v in obj.values():
            found.extend(walk_for_rows(v))
    return found


def parse_row(user, ts, cell):
    """Normalize one leaderboard row into a dict."""
    cell_unescaped = cell.encode().decode("unicode_escape") \
        if "\\u" in cell else cell
    if "Not submitted" in cell or "not submitted" in cell.lower():
        url = None
        status = "not_submitted"
    else:
        m = HREF_RE.search(cell_unescaped)
        url = m.group(1) if m else None
        status = "submitted" if url else "unknown"
    return {"username": user, "timestamp": ts, "status": status, "url": url}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True,
                    help="Path to the saved /config JSON file")
    ap.add_argument("--out", default="submissions.csv")
    args = ap.parse_args()

    rows = extract_from_json(args.config)
    # Deduplicate by (username, url), keep order
    seen = set()
    uniq = []
    for r in rows:
        key = (r["username"], r["url"])
        if key not in seen:
            seen.add(key)
            uniq.append(r)

    submitted = [r for r in uniq if r["status"] == "submitted"]
    print(f"Total rows found:     {len(uniq)}")
    print(f"With method links:    {len(submitted)}")
    print(f"Not submitted/blank:  {len(uniq) - len(submitted)}")
    print()
    print("=== Participants WITH method links ===")
    for r in submitted:
        print(f"  {r['username']:20s} {r['timestamp']:20s} {r['url']}")

    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["username", "timestamp",
                                          "status", "url"])
        w.writeheader()
        w.writerows(uniq)
    print(f"\nWrote {args.out}")
    # Also write just the URLs, one per line, for easy fetching
    urls_path = Path(args.out).with_suffix(".urls.txt")
    with open(urls_path, "w") as f:
        for r in submitted:
            if r["url"]:
                f.write(r["url"] + "\n")
    print(f"Wrote {urls_path} ({len(submitted)} URLs)")


if __name__ == "__main__":
    main()
