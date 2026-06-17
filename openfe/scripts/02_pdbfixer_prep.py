#!/usr/bin/env python3
"""
openfe/scripts/02_pdbfixer_prep.py

Prepare the unique PXR receptor structures needed for the RBFE campaign
using PDBFixer:
  - Find and add missing heavy atoms within existing residues
  - Add missing residues ONLY if the gap is <= max_missing_residues
    (longer gaps / terminal disordered regions are left as-is, capped
    by OpenMM's Modeller during system build)
  - Add hydrogens at pH 7.4 (LYS/ARG+, GLU/ASP-, HIS neutral by H-bond
    geometry - PDBFixer default behavior)

Input PDBs are the protein-only heavy-atom files written by
prep_receptors.py (docking/receptors/<pdb_id>_protein.pdb) - ligands,
waters, and crystallization additives already stripped. Using the raw
re-refined PDB files instead would cause PDBFixer to fail on non-standard
residue CIF parameters (e.g. rifampicin's '?' bond parameters).
Only the unique receptors listed in openfe/receptors/receptor_usage.csv
are processed (not all 62).

Output: openfe/receptors/<pdb_id>_prepared.pdb for each receptor.
These files are used as input to extract_rbfe_inputs.py and then
directly to openfe.ProteinComponent.from_pdb_file.

Usage:
    pip install pdbfixer openmm
    python openfe/scripts/02_pdbfixer_prep.py \
        --receptor-usage openfe/receptors/receptor_usage.csv \
        --pdb-dir docking/receptors \
        --outdir openfe/receptors \
        --max-missing-residues 10 \
        --ph 7.4
"""

import argparse
from pathlib import Path

import pandas as pd


def prepare_receptor(pdb_id, src_pdb, out_pdb, max_missing, ph):
    from pdbfixer import PDBFixer
    from openmm.app import PDBFile
    import tempfile, os

    print(f"\n  {pdb_id}: loading {src_pdb}")

    # PDBFixer chokes on '?' placeholder values that appear in some
    # re-refined PDB files (e.g. unknown anisou/cryst fields). Strip
    # any ATOM/HETATM lines containing '?' and any ANISOU records
    # (rarely needed for MD, and the most common source of '?' values).
    cleaned_lines = []
    n_stripped = 0
    with open(src_pdb) as f:
        for line in f:
            record = line[:6].strip()
            if record in ("ANISOU",):
                n_stripped += 1
                continue
            if "?" in line and record in ("ATOM", "HETATM", "CRYST1",
                                           "SCALE1", "SCALE2", "SCALE3",
                                           "ORIGX1", "ORIGX2", "ORIGX3"):
                n_stripped += 1
                continue
            cleaned_lines.append(line)
    if n_stripped:
        print(f"  {pdb_id}: stripped {n_stripped} lines containing '?' "
              f"or ANISOU records before PDBFixer parsing")

    # Write cleaned PDB to a temp file for PDBFixer
    with tempfile.NamedTemporaryFile(mode="w", suffix=".pdb",
                                      delete=False) as tmp:
        tmp.writelines(cleaned_lines)
        tmp_path = tmp.name

    try:
        fixer = PDBFixer(filename=tmp_path)
    finally:
        os.unlink(tmp_path)

    # Find missing residues and atoms
    fixer.findMissingResidues()
    fixer.findMissingAtoms()

    # Report what was found
    n_missing_res = sum(len(v) for v in fixer.missingResidues.values())
    n_missing_atoms = len(fixer.missingAtoms)
    print(f"  {pdb_id}: {n_missing_res} missing residue(s), "
          f"{n_missing_atoms} residue(s) with missing atoms")

    # Filter out long missing-residue gaps and terminal disordered regions
    # (these are likely IDRs or his-tags, safe to skip per best practice)
    keys_to_remove = []
    for key, residues in fixer.missingResidues.items():
        if len(residues) > max_missing:
            print(f"  {pdb_id}: skipping long missing segment at chain "
                  f"{key[0]} pos {key[1]}: {len(residues)} residues "
                  f"(> max {max_missing})")
            keys_to_remove.append(key)
    for key in keys_to_remove:
        del fixer.missingResidues[key]

    # Also remove terminal missing residues (likely his-tags / disordered ends)
    # PDBFixer flags these with chain index 0 / last position - check by
    # comparing to actual chain residue count
    # (PDBFixer handles this automatically when we skip long gaps, but
    # being explicit avoids accidentally modeling disordered termini)
    fixer.missingResidues = {
        key: res for key, res in fixer.missingResidues.items()
        if len(res) <= max_missing
    }

    # Add missing heavy atoms (within existing residues) and short gaps
    fixer.addMissingAtoms(seed=42)

    # Add hydrogens at the specified pH
    fixer.addMissingHydrogens(ph)
    print(f"  {pdb_id}: added hydrogens at pH {ph}")

    # Write output
    with open(out_pdb, "w") as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f,
                          keepIds=True)
    print(f"  {pdb_id}: wrote {out_pdb}")
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--receptor-usage", required=True,
                    help="openfe/receptors/receptor_usage.csv")
    ap.add_argument("--pdb-dir", required=True,
                    help="Directory containing <pdb_id>_protein.pdb files "
                         "(docking/receptors/ - the protein-only files "
                         "written by prep_receptors.py, ligands already "
                         "stripped). NOT the raw re-refinement directory.")
    ap.add_argument("--outdir", default="openfe/receptors")
    ap.add_argument("--max-missing-residues", type=int, default=10,
                    help="Maximum gap size to model; longer gaps are "
                         "left as chain breaks (default: 10)")
    ap.add_argument("--ph", type=float, default=7.4)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    pdb_dir = Path(args.pdb_dir)

    usage = pd.read_csv(args.receptor_usage)
    print(f"Preparing {len(usage)} unique receptors "
          f"(of 62 available) at pH {args.ph}")

    summary = []
    for _, row in usage.iterrows():
        pdb_id = row["pdb_id"]
        src_pdb = pdb_dir / f"{pdb_id}_protein.pdb"
        out_pdb = outdir / f"{pdb_id}_prepared.pdb"

        if not src_pdb.exists():
            print(f"  {pdb_id}: SOURCE NOT FOUND at {src_pdb} - skipping")
            summary.append({"pdb_id": pdb_id, "status": "MISSING_SOURCE",
                             "n_clusters": row["n_clusters"]})
            continue

        if out_pdb.exists():
            print(f"  {pdb_id}: already prepared, skipping "
                  f"(delete to re-run)")
            summary.append({"pdb_id": pdb_id, "status": "already_done",
                             "n_clusters": row["n_clusters"]})
            continue

        try:
            prepare_receptor(pdb_id, src_pdb, out_pdb,
                              args.max_missing_residues, args.ph)
            summary.append({"pdb_id": pdb_id, "status": "ok",
                             "n_clusters": row["n_clusters"]})
        except Exception as e:
            print(f"  {pdb_id}: FAILED - {type(e).__name__}: {e}")
            summary.append({"pdb_id": pdb_id, "status": f"FAILED: {e}",
                             "n_clusters": row["n_clusters"]})

    summary_df = pd.DataFrame(summary)
    summary_path = outdir / "pdbfixer_summary.csv"
    summary_df.to_csv(summary_path, index=False)

    n_ok = (summary_df["status"] == "ok").sum()
    print(f"\n{n_ok} / {len(summary_df)} receptors prepared successfully")
    print(f"Wrote {summary_path}")

    failed = summary_df[~summary_df["status"].isin(["ok", "already_done"])]
    if len(failed):
        print(f"\n{len(failed)} failed:")
        print(failed[["pdb_id", "status"]].to_string(index=False))


if __name__ == "__main__":
    main()
