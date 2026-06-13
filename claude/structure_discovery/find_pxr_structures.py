"""
PXR Challenge: locate and triage protein structure resources.

Run this LOCALLY (needs internet access to huggingface.co, github.com,
and files.rcsb.org).

This script checks THREE possible sources for PXR structures, since the
"68 re-refined PDB structures" mentioned in the challenge announcement
do not appear to be in the main pxr-challenge-train-test dataset repo:

  1. The 'structure' config of openadmet/pxr-challenge-train-test
     (78-compound structure TEST set - blinded ligands, but may
     reveal source PDB metadata / protein-only files)
  2. The OpenADMET/PXR-Challenge-Tutorial GitHub repo (may bundle
     re-refined PDBs and/or Boltz-2 reference structures)
  3. RCSB PDB fallback: a curated list of public human PXR LBD
     structures, downloaded directly if (1) and (2) come up empty.

Usage:
    pip install huggingface_hub requests pandas biopython
    python find_pxr_structures.py
"""

import os
import re
import json
import shutil
import zipfile
from pathlib import Path

import requests
import pandas as pd

OUT_DIR = Path("./pxr_structure_search")
OUT_DIR.mkdir(exist_ok=True)

report = {"hf_structure_config": None, "github_tutorial": None, "rcsb_fallback": None}

# ----------------------------------------------------------------------
# SOURCE 1: HF dataset 'structure' config
# ----------------------------------------------------------------------
print("=" * 70)
print("SOURCE 1: openadmet/pxr-challenge-train-test 'structure' config")
print("=" * 70)

from huggingface_hub import list_repo_files, hf_hub_download

REPO_ID = "openadmet/pxr-challenge-train-test"
try:
    all_files = list_repo_files(REPO_ID, repo_type="dataset")
    print(f"\nAll {len(all_files)} files in {REPO_ID}:")
    for f in all_files:
        print("  -", f)

    struct_files = [f for f in all_files if "structure" in f.lower()]
    print(f"\nFiles matching 'structure': {struct_files}")

    hf_struct_dir = OUT_DIR / "hf_structure_config"
    hf_struct_dir.mkdir(exist_ok=True)
    downloaded_struct = []
    for f in struct_files:
        local_path = hf_hub_download(
            repo_id=REPO_ID, repo_type="dataset", filename=f,
            local_dir=str(hf_struct_dir)
        )
        downloaded_struct.append(local_path)
        print("  downloaded ->", local_path)

    # if it's a CSV, peek at columns - look for PDB ID / structure refs
    for p in downloaded_struct:
        if p.endswith(".csv"):
            df = pd.read_csv(p)
            print(f"\n  {p} columns: {list(df.columns)}")
            print(f"  shape: {df.shape}")
            print(df.head(3).to_string())
            df.to_csv(OUT_DIR / f"hf_{Path(p).name}", index=False)

    report["hf_structure_config"] = {
        "files": struct_files,
        "downloaded": downloaded_struct,
    }
except Exception as e:
    print(f"  ERROR: {e}")
    report["hf_structure_config"] = {"error": str(e)}

# ----------------------------------------------------------------------
# SOURCE 2: GitHub tutorial repo
# ----------------------------------------------------------------------
print("\n" + "=" * 70)
print("SOURCE 2: OpenADMET/PXR-Challenge-Tutorial GitHub repo")
print("=" * 70)

GH_REPO = "OpenADMET/PXR-Challenge-Tutorial"
gh_dir = OUT_DIR / "github_tutorial"
gh_dir.mkdir(exist_ok=True)

try:
    # Get repo tree via GitHub API (default branch)
    api_url = f"https://api.github.com/repos/{GH_REPO}"
    resp = requests.get(api_url, timeout=30)
    resp.raise_for_status()
    default_branch = resp.json()["default_branch"]
    print(f"\nDefault branch: {default_branch}")

    tree_url = (f"https://api.github.com/repos/{GH_REPO}/git/trees/"
                f"{default_branch}?recursive=1")
    resp = requests.get(tree_url, timeout=30)
    resp.raise_for_status()
    tree = resp.json()["tree"]

    print(f"\nTotal files in repo: {len(tree)}")
    structure_like = [
        t["path"] for t in tree
        if re.search(r"\.(pdb|cif|ent|pdbqt)(\.gz)?$", t["path"], re.IGNORECASE)
        or re.search(r"(structure|boltz|pdb|crystal)", t["path"], re.IGNORECASE)
    ]
    print(f"\nFiles matching structure/pdb/boltz/crystal pattern "
          f"({len(structure_like)}):")
    for p in structure_like:
        print("  -", p)

    # Download matching files (small ones only - skip huge binaries)
    downloaded_gh = []
    for path in structure_like:
        raw_url = (f"https://raw.githubusercontent.com/{GH_REPO}/"
                   f"{default_branch}/{path}")
        try:
            r = requests.get(raw_url, timeout=30)
            if r.status_code == 200 and len(r.content) < 20_000_000:
                local_path = gh_dir / Path(path).name
                local_path.parent.mkdir(parents=True, exist_ok=True)
                with open(local_path, "wb") as f:
                    f.write(r.content)
                downloaded_gh.append(str(local_path))
                print(f"  downloaded {path} ({len(r.content)} bytes)")
            else:
                print(f"  skipped {path} (status {r.status_code}, "
                      f"size {len(r.content)})")
        except Exception as e:
            print(f"  failed {path}: {e}")

    report["github_tutorial"] = {
        "structure_like_files": structure_like,
        "downloaded": downloaded_gh,
    }
except Exception as e:
    print(f"  ERROR: {e}")
    report["github_tutorial"] = {"error": str(e)}

# ----------------------------------------------------------------------
# SOURCE 3: RCSB PDB fallback - public human PXR LBD structures
# ----------------------------------------------------------------------
print("\n" + "=" * 70)
print("SOURCE 3: RCSB PDB fallback (public human PXR LBD structures)")
print("=" * 70)

# A curated set of well-known human PXR (NR1I2) ligand-binding-domain
# structures from the PDB, spanning a range of bound ligand chemotypes.
# This list is NOT exhaustive - if sources 1/2 above yield the official
# 68-structure re-refined set, prefer that. This is a fallback to get
# started immediately.
PXR_PDB_IDS = [
    "1ILG",  # apo-like / SRC-1 peptide, early PXR LBD structure
    "1NRL",  # PXR with SR12813
    "1M13",  # PXR with hyperforin
    "1SKX",  # PXR with colupulone
    "2O9I",  # PXR with rifampicin-related ligand
    "1ILH",  # PXR with SRC-1 peptide
    "2QNV",  # PXR with T-0901317
    "3R8D",  # PXR with a synthetic ligand (mentioned in challenge text)
    "5X0Q",  # PXR with sulfasalazine-type ligand
    "6IBQ",  # PXR with an unspecified ligand
]

rcsb_dir = OUT_DIR / "rcsb_pdbs"
rcsb_dir.mkdir(exist_ok=True)

downloaded_rcsb = []
for pdb_id in PXR_PDB_IDS:
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    try:
        r = requests.get(url, timeout=30)
        if r.status_code == 200:
            local_path = rcsb_dir / f"{pdb_id}.pdb"
            with open(local_path, "wb") as f:
                f.write(r.content)
            downloaded_rcsb.append(str(local_path))
            print(f"  downloaded {pdb_id}.pdb ({len(r.content)} bytes)")
        else:
            print(f"  {pdb_id}: HTTP {r.status_code} (may not exist / "
                  f"check ID)")
    except Exception as e:
        print(f"  {pdb_id}: failed - {e}")

report["rcsb_fallback"] = {"attempted": PXR_PDB_IDS, "downloaded": downloaded_rcsb}

# ----------------------------------------------------------------------
# Write report
# ----------------------------------------------------------------------
with open(OUT_DIR / "search_report.json", "w") as f:
    json.dump(report, f, indent=2)

print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
print(f"Results written to {OUT_DIR}/")
print("Upload the following for further analysis:")
print(f"  - {OUT_DIR}/search_report.json")
print(f"  - any CSVs in {OUT_DIR}/ (HF structure config peek)")
print(f"  - {OUT_DIR}/rcsb_pdbs/*.pdb (fallback structures, if HF/GitHub empty)")
print(f"  - {OUT_DIR}/github_tutorial/* (if any structure files found)")

# Bundle everything for easy upload
shutil.make_archive(str(OUT_DIR.parent / "pxr_structure_search_bundle"), "zip",
                     root_dir=str(OUT_DIR))
print(f"\nBundled everything -> {OUT_DIR.parent / 'pxr_structure_search_bundle.zip'}")
