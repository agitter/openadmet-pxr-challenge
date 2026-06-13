"""
PXR Challenge: Download and triage the 68 re-refined PDB structures.

Run this LOCALLY (needs huggingface_hub + internet access to huggingface.co).
Outputs a CSV summary + a folder of extracted PDB files that can be
uploaded back for further analysis (pocket clustering, conformer selection).

Usage:
    pip install huggingface_hub rdkit biopython pandas
    python analyze_pxr_structures.py
"""

import os
import re
import json
import zipfile
import shutil
from pathlib import Path

import pandas as pd

# ----------------------------------------------------------------------
# STEP 1: Locate and download the structure training data package
# ----------------------------------------------------------------------
# The HF dataset repo is openadmet/pxr-challenge-train-test (same repo as
# the activity CSVs). The structure track training package (68 re-refined
# PDBs) should live alongside the other files. We list the repo contents
# first so this script is robust to filename changes.

from huggingface_hub import HfApi, hf_hub_download, list_repo_files

REPO_ID = "openadmet/pxr-challenge-train-test"
REPO_TYPE = "dataset"

OUT_DIR = Path("./pxr_structure_analysis")
PDB_DIR = OUT_DIR / "pdbs"
OUT_DIR.mkdir(exist_ok=True)
PDB_DIR.mkdir(exist_ok=True, parents=True)

print(f"Listing files in {REPO_ID} ...")
all_files = list_repo_files(REPO_ID, repo_type=REPO_TYPE)

# Look for anything that looks like a structure package: zip, tar, or
# a folder of PDB/CIF files, or anything with "structure" / "pdb" in the name
structure_files = [
    f for f in all_files
    if re.search(r"(structure|pdb|crystal|68)", f, re.IGNORECASE)
]

print(f"\nFound {len(all_files)} total files in repo.")
print("Files matching structure/pdb/crystal/68 pattern:")
for f in structure_files:
    print("  -", f)

if not structure_files:
    print("\nNo obvious structure package found by name. Printing ALL files")
    print("so you can identify the right one manually:")
    for f in all_files:
        print("  -", f)
    raise SystemExit(
        "\nEdit STRUCTURE_FILES below with the correct path(s) from the "
        "list above, then re-run."
    )

# Download each matching file
downloaded = []
for f in structure_files:
    print(f"\nDownloading {f} ...")
    local_path = hf_hub_download(
        repo_id=REPO_ID, repo_type=REPO_TYPE, filename=f,
        local_dir=str(OUT_DIR)
    )
    downloaded.append(local_path)
    print("  -> saved to", local_path)

# ----------------------------------------------------------------------
# STEP 2: Extract any archives, collect all .pdb / .cif / .pdb.gz files
# ----------------------------------------------------------------------
def extract_archive(path):
    p = Path(path)
    if p.suffix == ".zip":
        with zipfile.ZipFile(p) as zf:
            zf.extractall(PDB_DIR)
            print(f"  extracted zip -> {PDB_DIR}")
    elif p.suffixes[-2:] == [".tar", ".gz"] or p.suffix == ".tgz":
        import tarfile
        with tarfile.open(p) as tf:
            tf.extractall(PDB_DIR)
            print(f"  extracted tar.gz -> {PDB_DIR}")
    elif p.suffix in (".pdb", ".cif", ".ent"):
        shutil.copy(p, PDB_DIR / p.name)
    elif p.suffixes[-1:] == [".gz"]:
        import gzip
        out = PDB_DIR / p.stem
        with gzip.open(p, "rb") as fin, open(out, "wb") as fout:
            shutil.copyfileobj(fin, fout)

for f in downloaded:
    extract_archive(f)

pdb_files = sorted(PDB_DIR.rglob("*.pdb")) + sorted(PDB_DIR.rglob("*.cif")) \
    + sorted(PDB_DIR.rglob("*.ent"))
print(f"\nTotal structure files found after extraction: {len(pdb_files)}")
for pf in pdb_files[:10]:
    print("  -", pf)
if len(pdb_files) > 10:
    print(f"  ... and {len(pdb_files) - 10} more")

if not pdb_files:
    raise SystemExit("No PDB/CIF files found after extraction - inspect "
                      f"{OUT_DIR} manually.")

# ----------------------------------------------------------------------
# STEP 3: Parse each structure - extract resolution, chains, hetero
# residues (candidate bound ligands), and pocket-lining residue coords
# for later conformational clustering.
# ----------------------------------------------------------------------
from Bio.PDB import PDBParser, MMCIFParser
from Bio.PDB.Polypeptide import is_aa
import numpy as np
import warnings
from Bio import BiopythonWarning
warnings.simplefilter("ignore", BiopythonWarning)

# Standard amino acid 3-letter codes + common modified residues to exclude
# from "ligand" candidates
STANDARD_RESIDUES = set("""
ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL
HOH WAT NA CL MG ZN K CA SO4 PO4 GOL EDO DMS PEG ACT FMT MES TRS
""".split())

def parse_structure(path):
    path = Path(path)
    pdb_id = path.stem.upper()[:4]
    try:
        if path.suffix == ".cif":
            parser = MMCIFParser(QUIET=True)
        else:
            parser = PDBParser(QUIET=True)
        structure = parser.get_structure(pdb_id, str(path))
    except Exception as e:
        return {"file": str(path), "pdb_id": pdb_id, "error": str(e)}

    # Resolution: try header dict
    resolution = None
    try:
        if hasattr(structure, "header"):
            resolution = structure.header.get("resolution")
    except Exception:
        pass

    # Find hetero residues (candidate ligands) and chains
    hetero_residues = []
    chains = set()
    n_protein_residues = 0
    for model in structure:
        for chain in model:
            chains.add(chain.id)
            for residue in chain:
                resname = residue.get_resname().strip()
                het_flag = residue.id[0]
                if het_flag == " " or het_flag == "":
                    if is_aa(residue, standard=True):
                        n_protein_residues += 1
                elif het_flag.startswith("H_") or het_flag == "H":
                    if resname not in STANDARD_RESIDUES:
                        hetero_residues.append((chain.id, resname, residue.id))
                elif resname not in STANDARD_RESIDUES and resname != "HOH":
                    hetero_residues.append((chain.id, resname, residue.id))
        break  # only first model

    return {
        "file": str(path),
        "pdb_id": pdb_id,
        "resolution": resolution,
        "n_chains": len(chains),
        "chains": ",".join(sorted(chains)),
        "n_protein_residues": n_protein_residues,
        "n_hetero_residues": len(hetero_residues),
        "hetero_residue_names": ",".join(sorted(set(r[1] for r in hetero_residues))),
        "error": None,
    }


# Pocket definition: residues within 8A of ANY hetero atom that is NOT
# water/ion/crystallization additive (i.e. a putative ligand). We record
# CA coordinates of a fixed reference residue set for later clustering -
# but since we don't know the exact PXR numbering here, we instead just
# dump CA coordinates of all residues within 10A of the largest hetero
# group, for downstream alignment/clustering.

CRYO_ADDITIVES = {"GOL", "EDO", "DMS", "PEG", "ACT", "FMT", "MES", "TRS",
                  "SO4", "PO4", "NA", "CL", "MG", "ZN", "K", "CA", "HOH", "WAT"}

def get_pocket_ca_coords(path, cutoff=10.0):
    path = Path(path)
    pdb_id = path.stem.upper()[:4]
    try:
        if path.suffix == ".cif":
            parser = MMCIFParser(QUIET=True)
        else:
            parser = PDBParser(QUIET=True)
        structure = parser.get_structure(pdb_id, str(path))
    except Exception:
        return None, None

    model = next(iter(structure))

    # find largest non-additive hetero group -> treat as "ligand"
    ligand_atoms = []
    best_lig = None
    best_size = 0
    for chain in model:
        for residue in chain:
            resname = residue.get_resname().strip()
            het_flag = residue.id[0]
            is_het = het_flag != " " and het_flag != ""
            if is_het and resname not in CRYO_ADDITIVES and resname != "HOH":
                atoms = list(residue.get_atoms())
                if len(atoms) > best_size:
                    best_size = len(atoms)
                    best_lig = (chain.id, resname, residue.id)
                    ligand_atoms = atoms

    if not ligand_atoms:
        return None, None

    lig_coords = np.array([a.get_coord() for a in ligand_atoms])

    # pocket CA atoms within cutoff of any ligand atom
    pocket = []
    for chain in model:
        for residue in chain:
            if not is_aa(residue, standard=True):
                continue
            if "CA" not in residue:
                continue
            ca = residue["CA"].get_coord()
            dists = np.linalg.norm(lig_coords - ca, axis=1)
            if dists.min() <= cutoff:
                pocket.append((chain.id, residue.id[1], residue.get_resname(), ca))

    return best_lig, pocket


print("\nParsing structures ...")
records = []
pocket_data = {}
for pf in pdb_files:
    rec = parse_structure(pf)
    records.append(rec)
    lig, pocket = get_pocket_ca_coords(pf)
    rec["putative_ligand"] = f"{lig[0]}:{lig[1]}:{lig[2]}" if lig else None
    rec["n_pocket_residues"] = len(pocket) if pocket else 0
    pocket_data[rec["pdb_id"]] = {
        "file": str(pf),
        "ligand": lig,
        "pocket": [(c, num, name, coord.tolist()) for c, num, name, coord in (pocket or [])],
    }

df = pd.DataFrame(records)
df.to_csv(OUT_DIR / "structure_summary.csv", index=False)
with open(OUT_DIR / "pocket_data.json", "w") as f:
    json.dump(pocket_data, f, indent=2)

print("\n=== SUMMARY ===")
print(df[["pdb_id", "resolution", "n_chains", "n_protein_residues",
          "hetero_residue_names", "putative_ligand", "n_pocket_residues"]].to_string(index=False))

print(f"\nWrote summary -> {OUT_DIR / 'structure_summary.csv'}")
print(f"Wrote pocket coordinate data -> {OUT_DIR / 'pocket_data.json'}")
print(f"\nNow zip and upload the '{OUT_DIR}' folder (or at minimum "
      f"structure_summary.csv and pocket_data.json) for further analysis.")

# Also bundle PDB files for upload (optional - may be large)
shutil.make_archive(str(OUT_DIR / "pxr_structures_bundle"), "zip", root_dir=str(PDB_DIR))
print(f"Bundled PDB files -> {OUT_DIR / 'pxr_structures_bundle.zip'}")
