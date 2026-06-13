#!/usr/bin/env python3
"""
docking/scripts/prep_receptors.py

Prepare GNINA/AutoDock-ready receptor files from the 64 re-refined PXR
structures (external/pxr_xtal_re-refinement submodule).

For each structure:
  - Strip all HETATM records (ligand, waters, crystallization additives)
    leaving only the protein.
  - Save as a clean PDB (protein-only).
  - Convert to PDBQT using Meeko (preferred) or OpenBabel as fallback.
  - Record the docking box center (centroid of the original ligand's
    heavy atoms) and a default box size, written to receptors/boxes.csv.

Run this LOCALLY (once) - it does not need to run on CHTC. Output
(receptors/*.pdbqt + receptors/boxes.csv) gets staged to CHTC as input
for the docking jobs.

Usage:
    pip install biopython rdkit meeko
    python scripts/prep_receptors.py
"""

import csv
import warnings
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser, PDBIO, Select
from Bio.PDB.Polypeptide import is_aa
from Bio import BiopythonWarning

warnings.simplefilter("ignore", BiopythonWarning)

REPO_ROOT = Path("external/pxr_xtal_re-refinement/pxr_rerefined_structures")
OUT_DIR = Path("docking/receptors")
OUT_DIR.mkdir(parents=True, exist_ok=True)

EXCLUDE_RESNAMES = {
    "HOH", "WAT", "GOL", "EDO", "DMS", "PEG", "ACT", "FMT", "MES", "TRS",
    "SO4", "PO4", "NA", "CL", "MG", "ZN", "K", "CA", "1PE", "PG4", "BME",
    "IPA", "MPD", "PGE", "EPE", "CIT", "P6G", "OLA", "OLC", "OCT",
}

DEFAULT_BOX_SIZE = 25.0  # Angstrom, per dimension - generous for a large
                         # flexible pocket like PXR's


class ProteinOnly(Select):
    """Keep only standard amino acid residues (no waters, ligands, ions)."""
    def accept_residue(self, residue):
        return is_aa(residue, standard=True)

    def accept_atom(self, atom):
        # drop alternate conformers beyond the first (avoids duplicate-atom
        # issues in PDBQT conversion); keep atoms with altloc ' ' or 'A'
        altloc = atom.get_altloc()
        return altloc in (" ", "A")


def find_ligand_centroid(structure):
    """Find the largest non-excluded hetero residue and return the
    centroid of its heavy atoms, plus its identity."""
    model = next(iter(structure))
    best = None
    best_n = 0
    for chain in model:
        for residue in chain:
            resname = residue.get_resname().strip()
            het_flag = residue.id[0]
            is_het = het_flag != " " and het_flag != ""
            if is_het and resname not in EXCLUDE_RESNAMES:
                heavy_atoms = [a for a in residue.get_atoms() if a.element != "H"]
                if len(heavy_atoms) > best_n:
                    best_n = len(heavy_atoms)
                    best = (resname, heavy_atoms)
    if best is None:
        return None, None
    resname, atoms = best
    coords = np.array([a.get_coord() for a in atoms])
    centroid = coords.mean(axis=0)
    return resname, centroid


def main():
    pdb_ids = sorted(p.name for p in REPO_ROOT.iterdir() if p.is_dir())
    print(f"Found {len(pdb_ids)} structures")

    box_rows = []
    parser = PDBParser(QUIET=True)
    io = PDBIO()

    for pdb_id in pdb_ids:
        struct_dir = REPO_ROOT / pdb_id
        pdb_file = struct_dir / f"{pdb_id}.pdb"
        if not pdb_file.exists():
            print(f"  {pdb_id}: no .pdb file, skipping (use .cif manually if needed)")
            continue

        structure = parser.get_structure(pdb_id, str(pdb_file))

        # find ligand centroid for box placement BEFORE stripping
        lig_resname, centroid = find_ligand_centroid(structure)
        if centroid is None:
            print(f"  {pdb_id}: WARNING no ligand found for box centroid, "
                  f"using protein center of mass instead")
            model = next(iter(structure))
            coords = []
            for chain in model:
                for residue in chain:
                    if is_aa(residue, standard=True) and "CA" in residue:
                        coords.append(residue["CA"].get_coord())
            centroid = np.array(coords).mean(axis=0)
            lig_resname = None

        # write protein-only PDB
        protein_pdb = OUT_DIR / f"{pdb_id}_protein.pdb"
        io.set_structure(structure)
        io.save(str(protein_pdb), ProteinOnly())

        box_rows.append({
            "pdb_id": pdb_id,
            "ligand_resname": lig_resname,
            "center_x": round(float(centroid[0]), 3),
            "center_y": round(float(centroid[1]), 3),
            "center_z": round(float(centroid[2]), 3),
            "size_x": DEFAULT_BOX_SIZE,
            "size_y": DEFAULT_BOX_SIZE,
            "size_z": DEFAULT_BOX_SIZE,
        })
        print(f"  {pdb_id}: ligand={lig_resname}, "
              f"center=({centroid[0]:.1f},{centroid[1]:.1f},{centroid[2]:.1f}) "
              f"-> {protein_pdb}")

    with open(OUT_DIR / "boxes.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(box_rows[0].keys()))
        writer.writeheader()
        writer.writerows(box_rows)

    print(f"\nWrote {len(box_rows)} protein PDBs + box definitions to "
          f"{OUT_DIR}/boxes.csv")

    # -------------------------------------------------------------
    # PDBQT conversion
    # -------------------------------------------------------------
    print("\nConverting to PDBQT ...")
    converted = 0
    failed = []
    try:
        from meeko import MoleculePreparation, PDBQTWriterLegacy
        from rdkit import Chem
        use_meeko = True
    except ImportError:
        use_meeko = False
        print("  meeko not available, falling back to OpenBabel "
              "(pip install meeko for better results, or ensure "
              "`obabel` is on PATH)")

    import subprocess
    for row in box_rows:
        pdb_id = row["pdb_id"]
        protein_pdb = OUT_DIR / f"{pdb_id}_protein.pdb"
        out_pdbqt = OUT_DIR / f"{pdb_id}_protein.pdbqt"

        # Receptor PDBQT generation is most robustly done via
        # OpenBabel or AutoDockTools' prepare_receptor4.py (from
        # ADFR suite). Meeko is primarily for ligands. We use obabel
        # here for the receptor; ensure polar hydrogens are added.
        cmd = ["obabel", str(protein_pdb), "-O", str(out_pdbqt),
               "-xr", "--partialcharge", "gasteiger"]
        try:
            result = subprocess.run(cmd, capture_output=True, text=True,
                                     timeout=120)
            if result.returncode == 0 and out_pdbqt.exists():
                converted += 1
            else:
                failed.append((pdb_id, result.stderr[-300:]))
        except FileNotFoundError:
            failed.append((pdb_id, "obabel not found on PATH"))
        except Exception as e:
            failed.append((pdb_id, str(e)))

    print(f"\nConverted {converted}/{len(box_rows)} receptors to PDBQT")
    if failed:
        print("Failed conversions:")
        for pdb_id, err in failed:
            print(f"  {pdb_id}: {err}")
        print("\nIf obabel is unavailable, install it: "
              "conda install -c conda-forge openbabel, "
              "or apt-get install openbabel")


if __name__ == "__main__":
    main()
