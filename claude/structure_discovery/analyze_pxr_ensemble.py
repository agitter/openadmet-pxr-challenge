"""
PXR Challenge: analyze the 64 re-refined PXR crystal structures.

Run this LOCALLY from the root of your repo (where
external/pxr_xtal_re-refinement is checked out as a submodule).

For each structure, this script:
  - Parses the .pdb file (falls back to .cif if .pdb absent, e.g. 9fzg/9fzh)
  - Identifies the bound ligand (largest non-water, non-crystallization-
    additive HETATM residue)
  - Extracts the ligand's SMILES via RDKit (from PDB connectivity - may
    need manual correction for structures where bond orders are ambiguous)
  - Extracts pocket residues (protein residues with any atom within 8A
    of any ligand atom) - residue identity + CA coordinates
  - Computes ligand heavy-atom count, formula, and basic descriptors

Output: pxr_structure_inventory.csv (one row per structure) and
        pxr_pocket_residues.json (pocket residue lists + CA coords per
        structure, for later pocket-similarity clustering)

Usage:
    pip install rdkit biopython pandas
    python analyze_pxr_ensemble.py
"""

import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from Bio.PDB import PDBParser, MMCIFParser
from Bio.PDB.Polypeptide import is_aa
from Bio import BiopythonWarning

warnings.simplefilter("ignore", BiopythonWarning)

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

REPO_ROOT = Path("external/pxr_xtal_re-refinement/pxr_rerefined_structures")

# Known crystallization additives / ions / waters to exclude when looking
# for the "real" bound ligand
EXCLUDE_RESNAMES = {
    "HOH", "WAT", "GOL", "EDO", "DMS", "PEG", "ACT", "FMT", "MES", "TRS",
    "SO4", "PO4", "NA", "CL", "MG", "ZN", "K", "CA", "1PE", "PG4", "BME",
    "IPA", "MPD", "PGE", "EPE", "CIT", "P6G", "OLA", "OLC", "OCT",
}

STANDARD_AA = set("""
ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL
""".split())


def get_parser(path: Path):
    if path.suffix.lower() == ".cif":
        return MMCIFParser(QUIET=True)
    return PDBParser(QUIET=True)


def find_ligand_residues(structure):
    """Return list of (chain_id, resname, res_id, residue) for candidate
    ligand hetero residues, sorted by heavy-atom count descending."""
    candidates = []
    model = next(iter(structure))
    for chain in model:
        for residue in chain:
            resname = residue.get_resname().strip()
            het_flag = residue.id[0]
            is_het = het_flag != " " and het_flag != ""
            if is_het and resname not in EXCLUDE_RESNAMES:
                n_heavy = sum(1 for a in residue.get_atoms()
                               if a.element != "H")
                if n_heavy >= 5:  # filter out tiny fragments/ions we missed
                    candidates.append((chain.id, resname, residue.id, residue, n_heavy))
    candidates.sort(key=lambda x: -x[4])
    return candidates


def residue_to_mol(residue, pdb_id, resname):
    """Attempt to build an RDKit mol from a Biopython residue via a
    temporary PDB block. Returns (mol, smiles) or (None, None)."""
    from io import StringIO
    from Bio.PDB import PDBIO, Select

    class ResSelect(Select):
        def accept_residue(self, res):
            return res is residue

    pdbio = PDBIO()
    pdbio.set_structure(residue.get_parent().get_parent().get_parent())
    sio = StringIO()
    try:
        pdbio.save(sio, ResSelect())
    except Exception:
        return None, None
    pdb_block = sio.getvalue()
    if not pdb_block.strip():
        return None, None

    mol = Chem.MolFromPDBBlock(pdb_block, sanitize=False, removeHs=False)
    if mol is None:
        return None, None
    try:
        Chem.SanitizeMol(mol, sanitizeOps=Chem.SANITIZE_ALL ^ Chem.SANITIZE_ADJUSTHS)
    except Exception:
        pass
    try:
        smiles = Chem.MolToSmiles(mol)
    except Exception:
        smiles = None
    return mol, smiles


def get_pocket_residues(structure, ligand_residue, cutoff=8.0):
    """Return list of (chain_id, resnum, resname, ca_coord) for protein
    residues with any atom within cutoff Angstrom of any ligand atom."""
    lig_coords = np.array([a.get_coord() for a in ligand_residue.get_atoms()
                            if a.element != "H"])
    if len(lig_coords) == 0:
        return []

    model = next(iter(structure))
    pocket = []
    for chain in model:
        for residue in chain:
            if not is_aa(residue, standard=True):
                continue
            min_dist = np.inf
            for atom in residue.get_atoms():
                if atom.element == "H":
                    continue
                d = np.linalg.norm(lig_coords - atom.get_coord(), axis=1).min()
                if d < min_dist:
                    min_dist = d
            if min_dist <= cutoff:
                ca_coord = (residue["CA"].get_coord().tolist()
                            if "CA" in residue else None)
                pocket.append({
                    "chain": chain.id,
                    "resnum": residue.id[1],
                    "resname": residue.get_resname(),
                    "min_dist_to_ligand": round(float(min_dist), 2),
                    "ca_coord": ca_coord,
                })
    pocket.sort(key=lambda x: x["min_dist_to_ligand"])
    return pocket


def get_resolution(structure, path):
    """Try to extract resolution from header (PDB) - CIF parser in
    Biopython often doesn't populate this, so fall back to grepping
    the file for common resolution tags."""
    res = None
    try:
        res = structure.header.get("resolution")
    except Exception:
        pass
    if res is None:
        # crude grep fallback for mmCIF _refine.ls_d_res_high or
        # PDB REMARK 2 RESOLUTION
        try:
            text = path.read_text(errors="ignore")
            for line in text.splitlines():
                if line.startswith("REMARK   2 RESOLUTION"):
                    parts = line.split()
                    for p in parts:
                        try:
                            res = float(p)
                            break
                        except ValueError:
                            continue
                    if res:
                        break
                if "_refine.ls_d_res_high" in line:
                    parts = line.split()
                    if len(parts) > 1:
                        try:
                            res = float(parts[-1])
                        except ValueError:
                            pass
                    break
        except Exception:
            pass
    return res


def main():
    pdb_ids = sorted(p.name for p in REPO_ROOT.iterdir() if p.is_dir())
    print(f"Found {len(pdb_ids)} structure directories")

    records = []
    pocket_data = {}

    for pdb_id in pdb_ids:
        struct_dir = REPO_ROOT / pdb_id
        pdb_file = struct_dir / f"{pdb_id}.pdb"
        cif_file = struct_dir / f"{pdb_id}.cif"

        path = pdb_file if pdb_file.exists() else cif_file
        if not path.exists():
            print(f"  {pdb_id}: NO STRUCTURE FILE FOUND, skipping")
            continue

        parser = get_parser(path)
        try:
            structure = parser.get_structure(pdb_id, str(path))
        except Exception as e:
            print(f"  {pdb_id}: PARSE ERROR - {e}")
            records.append({"pdb_id": pdb_id, "error": str(e)})
            continue

        resolution = get_resolution(structure, path)

        candidates = find_ligand_residues(structure)
        if not candidates:
            print(f"  {pdb_id}: no ligand candidates found")
            records.append({
                "pdb_id": pdb_id, "file_used": path.name,
                "resolution": resolution,
                "ligand_resname": None, "ligand_smiles": None,
                "ligand_n_heavy": None, "n_pocket_residues": 0,
                "error": "no_ligand_found",
            })
            continue

        # take the largest candidate as the primary ligand
        chain_id, resname, res_id, residue, n_heavy = candidates[0]
        mol, smiles = residue_to_mol(residue, pdb_id, resname)

        mw = formula = None
        if mol is not None:
            try:
                mw = round(Descriptors.MolWt(mol), 1)
                formula = rdMolDescriptors.CalcMolFormula(mol)
            except Exception:
                pass

        pocket = get_pocket_residues(structure, residue, cutoff=8.0)

        other_ligands = [f"{r[1]}({r[4]}ha)" for r in candidates[1:]]

        records.append({
            "pdb_id": pdb_id,
            "file_used": path.name,
            "resolution": resolution,
            "ligand_chain": chain_id,
            "ligand_resname": resname,
            "ligand_resid": str(res_id),
            "ligand_n_heavy": n_heavy,
            "ligand_smiles": smiles,
            "ligand_mw": mw,
            "ligand_formula": formula,
            "n_other_ligand_candidates": len(candidates) - 1,
            "other_ligand_resnames": ",".join(other_ligands) if other_ligands else None,
            "n_pocket_residues": len(pocket),
            "error": None,
        })

        pocket_data[pdb_id] = {
            "ligand_resname": resname,
            "ligand_smiles": smiles,
            "pocket_residues": pocket,
        }

        print(f"  {pdb_id}: ligand={resname} ({n_heavy} heavy atoms, "
              f"formula={formula}), pocket_residues={len(pocket)}, "
              f"resolution={resolution}")

    df = pd.DataFrame(records)
    df.to_csv("pxr_structure_inventory.csv", index=False)
    with open("pxr_pocket_residues.json", "w") as f:
        json.dump(pocket_data, f, indent=2)

    print(f"\nWrote pxr_structure_inventory.csv ({len(df)} rows)")
    print(f"Wrote pxr_pocket_residues.json")

    n_with_smiles = df["ligand_smiles"].notna().sum()
    print(f"\nStructures with successfully extracted ligand SMILES: "
          f"{n_with_smiles} / {len(df)}")
    print("\nNOTE: SMILES extracted from PDB connectivity via RDKit can "
          "have WRONG BOND ORDERS / protonation states (PDB format does "
          "not reliably encode these). Treat 'ligand_smiles' as a rough "
          "identity hint for similarity matching, not a validated "
          "structure. Cross-check important ligands manually against "
          "PDB ligand IDs (ligand_resname) at "
          "https://www.rcsb.org/ligand/<resname>")


if __name__ == "__main__":
    main()
