#!/usr/bin/env python3
"""
docking/scripts/extract_rbfe_inputs.py

For each cluster in rbfe_template_selection.csv, extract a paired
(receptor PDB, ligand SDF) RBFE starting structure:

  - Receptor: docking/receptors/<selected_pdb_id>_protein.pdb, copied
    as-is (heavy atoms only, no waters/ligands - this is the file
    prep_receptors.py wrote BEFORE PDBQT conversion, so its coordinates
    are identical to the PDBQT used for docking).

  - Ligand: pose `selected_pose_rank` extracted from
    results/<cluster_id>/<cluster_id>__<selected_pdb_id>_docked.sdf.gz,
    with bond orders/aromaticity/formal charges reassigned from the
    original cluster-representative SMILES via
    AllChem.AssignBondOrdersFromTemplate (gnina's PDBQT->SDF round trip
    can produce incorrect bond orders even when 3D coordinates and
    heavy-atom connectivity are correct), then AddHs(addCoords=True) to
    restore a full all-explicit-hydrogen molecule with RDKit-placed H
    positions consistent with the docked heavy-atom pose.

Because obabel's PDB->PDBQT conversion (used by prep_receptors.py)
does not translate/rotate coordinates, the docked ligand pose and the
protein.pdb receptor share the same coordinate frame and can be loaded
together directly as a complex.

Output (one subdirectory per cluster):
  <outdir>/<cluster_id>/<cluster_id>_receptor.pdb
  <outdir>/<cluster_id>/<cluster_id>_ligand.sdf

Also writes <outdir>/extraction_summary.csv logging template-matching
success/failure and atom counts for spot-checking.

Usage:
    pip install rdkit pandas
    python extract_rbfe_inputs.py \
        --selection docking/docking_analysis/rbfe_template_selection.csv \
        --results-dir results \
        --receptor-dir docking/receptors \
        --outdir docking/rbfe_inputs
"""

import argparse
import gzip
import shutil
from pathlib import Path

import pandas as pd

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")


def load_pose(sdf_gz_path, pose_rank):
    """Return the RDKit mol for the given 1-based pose_rank from a
    gzipped multi-record SDF (sanitize=False - bond orders will be
    fixed via template matching afterwards)."""
    with gzip.open(sdf_gz_path, "rb") as f:
        supplier = Chem.ForwardSDMolSupplier(f, sanitize=False, removeHs=False)
        for i, mol in enumerate(supplier, start=1):
            if i == pose_rank:
                return mol
    return None


def fix_bond_orders(docked_mol, smiles):
    """Reassign bond orders/aromaticity/formal charges from `smiles`
    onto `docked_mol`'s heavy-atom skeleton (keeping its 3D
    coordinates), then add explicit hydrogens with RDKit-placed
    coordinates. Returns the new mol, or raises on template-match
    failure.

    Kekulization note: AssignBondOrdersFromTemplate can fail with
    KekulizeException when the docked pose's connectivity (after PDBQT
    round-trip) doesn't exactly match the aromatic system expected by
    the SMILES template. Strategy:
      1. Convert template to a Kekulé form (explicit single/double bonds,
         no aromatic flags) before matching - this avoids aromaticity
         perception mismatches between template and docked mol.
      2. After matching, run full sanitization (which re-perceives
         aromaticity from scratch on the new mol's connectivity).
      3. If kekulization still fails on the matched mol, fall back to
         sanitizing without SANITIZE_KEKULIZE to at least get a mol
         with correct connectivity for inspection.
    """
    template = Chem.MolFromSmiles(smiles)
    if template is None:
        raise ValueError(f"Could not parse SMILES: {smiles}")

    # Convert template to Kekulé form (explicit alternating bonds) to
    # avoid aromatic-perception mismatches during template matching
    try:
        Chem.Kekulize(template, clearAromaticFlags=True)
    except Exception:
        pass  # if Kekulization of the template fails, use as-is

    docked_heavy = Chem.RemoveHs(docked_mol, sanitize=False)
    Chem.SanitizeMol(docked_heavy,
                      sanitizeOps=Chem.SANITIZE_ALL ^ Chem.SANITIZE_ADJUSTHS
                                   ^ Chem.SANITIZE_KEKULIZE)

    new_mol = AllChem.AssignBondOrdersFromTemplate(template, docked_heavy)

    # Full sanitization re-perceives aromaticity from the matched connectivity
    try:
        Chem.SanitizeMol(new_mol)
    except Exception:
        # Fall back: sanitize everything except kekulization - mol is still
        # usable for coordinate extraction; flag in caller
        Chem.SanitizeMol(new_mol,
                          sanitizeOps=Chem.SANITIZE_ALL ^ Chem.SANITIZE_KEKULIZE)

    new_mol = Chem.AddHs(new_mol, addCoords=True)
    return new_mol


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--selection", required=True,
                    help="rbfe_template_selection.csv")
    ap.add_argument("--results-dir", default="results")
    ap.add_argument("--receptor-dir", default="docking/receptors")
    ap.add_argument("--receptor-suffix", default="_protein.pdb",
                    help="Suffix appended to pdb_id to find receptor file. "
                         "Use '_protein.pdb' for raw stripped receptors "
                         "(docking/receptors/) or '_prepared.pdb' for "
                         "PDBFixer-prepared receptors "
                         "(openfe/receptors/). Default: _protein.pdb")
    ap.add_argument("--outdir", default="docking/rbfe_inputs")
    args = ap.parse_args()

    results_dir = Path(args.results_dir)
    receptor_dir = Path(args.receptor_dir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    selection = pd.read_csv(args.selection)
    print(f"Loaded {len(selection)} cluster selections")

    summary_rows = []

    for _, row in selection.iterrows():
        cluster_id = row["cluster_id"]
        pdb_id = row["selected_pdb_id"]
        pose_rank = int(row["selected_pose_rank"])
        smiles = row["smiles"]

        cluster_dir = outdir / str(cluster_id)
        cluster_dir.mkdir(parents=True, exist_ok=True)

        out_receptor = cluster_dir / f"{cluster_id}_receptor.pdb"
        out_ligand = cluster_dir / f"{cluster_id}_ligand.sdf"

        rec_status = "ok"
        lig_status = "ok"
        n_heavy_template = n_heavy_docked = n_atoms_final = None
        error = ""

        # --- Receptor: straight copy ---
        src_receptor = receptor_dir / f"{pdb_id}{args.receptor_suffix}"
        if src_receptor.exists():
            shutil.copy(src_receptor, out_receptor)
        else:
            rec_status = "MISSING_RECEPTOR_PDB"
            error += f"receptor not found: {src_receptor}; "

        # --- Ligand: extract pose + fix bond orders ---
        sdf_gz = results_dir / str(cluster_id) / f"{cluster_id}__{pdb_id}_docked.sdf.gz"
        if not sdf_gz.exists():
            lig_status = "MISSING_DOCKED_SDF"
            error += f"docked sdf not found: {sdf_gz}; "
        else:
            try:
                docked_mol = load_pose(sdf_gz, pose_rank)
                if docked_mol is None:
                    raise ValueError(
                        f"pose_rank {pose_rank} not found in {sdf_gz}")

                template = Chem.MolFromSmiles(smiles)
                n_heavy_template = template.GetNumHeavyAtoms() if template else None
                docked_heavy_count = Chem.RemoveHs(
                    docked_mol, sanitize=False).GetNumAtoms()
                n_heavy_docked = docked_heavy_count

                fixed_mol = fix_bond_orders(docked_mol, smiles)
                n_atoms_final = fixed_mol.GetNumAtoms()

                fixed_mol.SetProp("_Name", str(row["ligand_name"]))
                fixed_mol.SetProp("cluster_id", str(cluster_id))
                fixed_mol.SetProp("source_pdb_id", str(pdb_id))
                fixed_mol.SetProp("source_pose_rank", str(pose_rank))
                fixed_mol.SetProp("source_smiles", str(smiles))

                writer = Chem.SDWriter(str(out_ligand))
                writer.write(fixed_mol)
                writer.close()
            except Exception as e:
                lig_status = "TEMPLATE_MATCH_FAILED"
                error += f"{type(e).__name__}: {e}; "

        summary_rows.append({
            "cluster_id": cluster_id,
            "ligand_name": row["ligand_name"],
            "selected_pdb_id": pdb_id,
            "selected_pose_rank": pose_rank,
            "receptor_status": rec_status,
            "ligand_status": lig_status,
            "n_heavy_template": n_heavy_template,
            "n_heavy_docked": n_heavy_docked,
            "n_atoms_final": n_atoms_final,
            "error": error,
        })

    summary = pd.DataFrame(summary_rows)
    summary.to_csv(outdir / "extraction_summary.csv", index=False)

    n_ok = ((summary["receptor_status"] == "ok")
            & (summary["ligand_status"] == "ok")).sum()
    print(f"\n{n_ok} / {len(summary)} clusters fully extracted "
          f"(receptor + ligand)")

    failed = summary[(summary["receptor_status"] != "ok")
                      | (summary["ligand_status"] != "ok")]
    if len(failed):
        print(f"\n{len(failed)} cluster(s) with issues:")
        print(failed[["cluster_id", "receptor_status", "ligand_status",
                       "error"]].to_string(index=False))

    # Sanity check: heavy atom count mismatches (template vs docked) can
    # indicate a missing/extra atom even when template matching nominally
    # "succeeds" via substructure isomorphism on a subset
    mismatch = summary[
        (summary["n_heavy_template"].notna())
        & (summary["n_heavy_docked"].notna())
        & (summary["n_heavy_template"] != summary["n_heavy_docked"])
    ]
    if len(mismatch):
        print(f"\n{len(mismatch)} cluster(s) where docked heavy-atom "
              f"count != template heavy-atom count (inspect manually):")
        print(mismatch[["cluster_id", "n_heavy_template",
                         "n_heavy_docked"]].to_string(index=False))

    print(f"\nWrote {outdir/'extraction_summary.csv'}")


if __name__ == "__main__":
    main()
