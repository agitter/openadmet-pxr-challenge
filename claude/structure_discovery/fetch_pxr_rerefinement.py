"""
PXR Challenge: fetch the OpenADMET PXR crystal re-refinement repo.

Run this LOCALLY (needs internet access to github.com / raw.githubusercontent.com
and api.github.com).

This repo (https://github.com/OpenADMET/pxr_xtal_re-refinement) is the source
referenced in the OpenADMET blog post "Pregnane X Receptor PDB Structure
Rerefinement" - it should contain the 66 re-refined PXR-ligand PDB structures
(or the scripts + PDB ID list needed to regenerate them) used as the basis
for the structure-track training data.

Usage:
    pip install requests pandas
    python fetch_pxr_rerefinement.py
"""

import json
import re
from pathlib import Path

import requests

OUT_DIR = Path("./pxr_rerefinement")
OUT_DIR.mkdir(exist_ok=True)

GH_REPO = "OpenADMET/pxr_xtal_re-refinement"

print("=" * 70)
print(f"Inspecting {GH_REPO}")
print("=" * 70)

# 1. Get repo metadata (default branch, description, size)
api_url = f"https://api.github.com/repos/{GH_REPO}"
resp = requests.get(api_url, timeout=30)
if resp.status_code != 200:
    print(f"ERROR fetching repo metadata: HTTP {resp.status_code}")
    print(resp.text[:1000])
    raise SystemExit("Could not access repo - check name/spelling or "
                      "whether it's private.")

meta = resp.json()
default_branch = meta["default_branch"]
print(f"\nDescription: {meta.get('description')}")
print(f"Default branch: {default_branch}")
print(f"Size (KB): {meta.get('size')}")
print(f"Stars: {meta.get('stargazers_count')}")

with open(OUT_DIR / "repo_metadata.json", "w") as f:
    json.dump(meta, f, indent=2)

# 2. Get full file tree (recursive)
tree_url = (f"https://api.github.com/repos/{GH_REPO}/git/trees/"
            f"{default_branch}?recursive=1")
resp = requests.get(tree_url, timeout=60)
resp.raise_for_status()
tree = resp.json()["tree"]

print(f"\nTotal files/dirs in repo: {len(tree)}")

# Categorize files
pdb_like = []
cif_like = []
restraint_like = []
scripts = []
readmes = []
csv_json = []
other = []

for t in tree:
    if t["type"] != "blob":
        continue
    path = t["path"]
    lower = path.lower()
    if re.search(r"\.(pdb|ent)(\.gz)?$", lower):
        pdb_like.append(t)
    elif re.search(r"\.(cif|mmcif)(\.gz)?$", lower):
        cif_like.append(t)
    elif re.search(r"(restraint|cif\.params|\.params|eff$|\.eff)", lower):
        restraint_like.append(t)
    elif re.search(r"\.(py|sh|ipynb|smk)$", lower):
        scripts.append(t)
    elif re.search(r"readme", lower):
        readmes.append(t)
    elif re.search(r"\.(csv|json|tsv|yaml|yml)$", lower):
        csv_json.append(t)
    else:
        other.append(t)

print(f"\nPDB-like files: {len(pdb_like)}")
for t in pdb_like[:20]:
    print("  -", t["path"], f"({t.get('size', '?')} bytes)")
if len(pdb_like) > 20:
    print(f"  ... and {len(pdb_like)-20} more")

print(f"\nCIF-like files: {len(cif_like)}")
for t in cif_like[:20]:
    print("  -", t["path"], f"({t.get('size', '?')} bytes)")
if len(cif_like) > 20:
    print(f"  ... and {len(cif_like)-20} more")

print(f"\nRestraint-like files: {len(restraint_like)}")
for t in restraint_like[:20]:
    print("  -", t["path"])

print(f"\nScripts: {len(scripts)}")
for t in scripts[:30]:
    print("  -", t["path"])

print(f"\nCSV/JSON/YAML (likely PDB ID lists / metadata): {len(csv_json)}")
for t in csv_json[:30]:
    print("  -", t["path"])

print(f"\nREADMEs: {len(readmes)}")
for t in readmes:
    print("  -", t["path"])

print(f"\nOther files ({len(other)}):")
for t in other[:40]:
    print("  -", t["path"])
if len(other) > 40:
    print(f"  ... and {len(other)-40} more")

with open(OUT_DIR / "file_tree.json", "w") as f:
    json.dump(tree, f, indent=2)

# 3. Download README(s) and any small CSV/JSON/YAML metadata files
# (these likely contain the PDB ID list and/or pointers to where the
# actual structure files live, e.g. a release or LFS/zenodo link)
download_dir = OUT_DIR / "downloaded"
download_dir.mkdir(exist_ok=True)

to_download = readmes + [t for t in csv_json if t.get("size", 0) < 2_000_000]
# also grab top-level scripts that might define the PDB ID list
to_download += [t for t in scripts if t.get("size", 0) < 200_000
                 and "/" not in t["path"]]

print(f"\nDownloading {len(to_download)} small metadata/readme/script files ...")
for t in to_download:
    path = t["path"]
    raw_url = f"https://raw.githubusercontent.com/{GH_REPO}/{default_branch}/{path}"
    try:
        r = requests.get(raw_url, timeout=30)
        if r.status_code == 200:
            local_path = download_dir / path.replace("/", "__")
            with open(local_path, "wb") as f:
                f.write(r.content)
            print(f"  downloaded {path}")
        else:
            print(f"  failed {path}: HTTP {r.status_code}")
    except Exception as e:
        print(f"  failed {path}: {e}")

# 4. Check for GitHub Releases (large structure bundles often live here,
# not in the git tree directly)
print("\n" + "=" * 70)
print("Checking GitHub Releases for structure bundle downloads")
print("=" * 70)
releases_url = f"https://api.github.com/repos/{GH_REPO}/releases"
resp = requests.get(releases_url, timeout=30)
if resp.status_code == 200:
    releases = resp.json()
    print(f"\nFound {len(releases)} release(s)")
    for rel in releases:
        print(f"\n  Release: {rel.get('tag_name')} - {rel.get('name')}")
        for asset in rel.get("assets", []):
            print(f"    asset: {asset['name']} "
                  f"({asset['size']} bytes) -> {asset['browser_download_url']}")
    with open(OUT_DIR / "releases.json", "w") as f:
        json.dump(releases, f, indent=2)
else:
    print(f"  No releases endpoint accessible (HTTP {resp.status_code})")

print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
print(f"Inspect {OUT_DIR}/file_tree.json for the full repo listing,")
print(f"and {OUT_DIR}/downloaded/ for README/metadata contents.")
print("Upload these (or at least file_tree.json + the downloaded/ folder + "
      "releases.json) for further analysis.")

import shutil
shutil.make_archive(str(OUT_DIR.parent / "pxr_rerefinement_bundle"), "zip",
                     root_dir=str(OUT_DIR))
print(f"\nBundled -> {OUT_DIR.parent / 'pxr_rerefinement_bundle.zip'}")
