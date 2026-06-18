#!/usr/bin/env python3
"""
openfe/scripts/03_extract_all_rbfe_inputs.py

Extract RBFE starting structures (receptor PDB + ligand SDF) for every
node in the RBFE perturbation networks:
  - All 513 test compounds (T### IDs from extended docking campaign)
  - All 89 unique training anchors (A#### IDs) - every training compound
    that is the nearest neighbor of at least one test compound. Some
    clusters have multiple anchors (up to 4); all are included to minimize
    path length from any test compound to a known-pEC50 anchor node.

Each compound is assigned its cluster's selected receptor (from
rbfe_template_selection.csv) so that all compounds in a cluster share
the same receptor conformation.

The receptor used is the PDBFixer-prepared version (<pdb_id>_prepared.pdb
from openfe/receptors/).

The ligand pose is extracted from the extended docking campaign results
using the cluster's selected pdb_id and pose_rank. Bond orders are
reassigned from SMILES via AllChem.AssignBondOrdersFromTemplate with
Kekulé-form template matching.

Output structure:
  openfe/rbfe_inputs/<cluster_id>/
      T###_receptor.pdb, T###_ligand.sdf     (test compounds)
      A####_receptor.pdb, A####_ligand.sdf   (training anchors)
  openfe/rbfe_inputs/extraction_summary_extended.csv
  openfe/rbfe_inputs/per_cluster_extraction_counts.csv

Usage:
    pip install rdkit pandas
    python openfe/scripts/03_extract_all_rbfe_inputs.py \
        --selection docking/docking_analysis/rbfe_template_selection.csv \
        --test-full claude/outputs/test_full_with_clusters_and_anchors.csv \
        --train data/pxr-challenge_TRAIN.csv \
        --results-dir docking/results \
        --receptor-dir openfe/receptors \
        --outdir openfe/rbfe_inputs
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
    """Return RDKit mol for 1-based pose_rank from a gzipped SDF."""
    with gzip.open(sdf_gz_path, "rb") as f:
        supplier = Chem.ForwardSDMolSupplier(f, sanitize=False,
                                              removeHs=False)
        for i, mol in enumerate(supplier, start=1):
            if i == pose_rank:
                return mol
    return None


def fix_bond_orders(docked_mol, smiles):
    """Reassign bond orders from SMILES onto docked pose coordinates.
    Uses Kekulé-form template to avoid aromaticity mismatches."""
    template = Chem.MolFromSmiles(smiles)
    if template is None:
        raise ValueError(f"Could not parse SMILES: {smiles}")
    try:
        Chem.Kekulize(template, clearAromaticFlags=True)
    except Exception:
        pass

    docked_heavy = Chem.RemoveHs(docked_mol, sanitize=False)
    Chem.SanitizeMol(docked_heavy,
                      sanitizeOps=Chem.SANITIZE_ALL
                                   ^ Chem.SANITIZE_ADJUSTHS
                                   ^ Chem.SANITIZE_KEKULIZE)

    new_mol = AllChem.AssignBondOrdersFromTemplate(template, docked_heavy)
    try:
        Chem.SanitizeMol(new_mol)
    except Exception:
        Chem.SanitizeMol(new_mol,
                          sanitizeOps=Chem.SANITIZE_ALL
                                       ^ Chem.SANITIZE_KEKULIZE)

    return Chem.AddHs(new_mol, addCoords=True)


def extract_one(ligand_id, ligand_name, smiles, pdb_id, pose_rank,
                 cluster_id, results_dir, receptor_dir, cluster_dir,
                 pec50=None):
    """Extract receptor PDB + ligand SDF for one compound."""
    out_receptor = cluster_dir / f"{ligand_id}_receptor.pdb"
    out_ligand = cluster_dir / f"{ligand_id}_ligand.sdf"

    rec_status = lig_status = "ok"
    n_heavy_template = n_heavy_docked = n_atoms_final = None
    error = ""

    src_receptor = receptor_dir / f"{pdb_id}_prepared.pdb"
    if src_receptor.exists():
        shutil.copy(src_receptor, out_receptor)
    else:
        rec_status = "MISSING_RECEPTOR_PDB"
        error += f"receptor not found: {src_receptor}; "

    sdf_gz = (results_dir / ligand_id
               / f"{ligand_id}__{pdb_id}_docked.sdf.gz")
    if not sdf_gz.exists():
        lig_status = "MISSING_DOCKED_SDF"
        error += f"docked sdf not found: {sdf_gz}; "
    else:
        try:
            docked_mol = load_pose(sdf_gz, pose_rank)
            if docked_mol is None:
                raise ValueError(f"pose_rank {pose_rank} not in {sdf_gz}")

            template = Chem.MolFromSmiles(smiles)
            n_heavy_template = (template.GetNumHeavyAtoms()
                                 if template else None)
            n_heavy_docked = Chem.RemoveHs(
                docked_mol, sanitize=False).GetNumAtoms()

            fixed_mol = fix_bond_orders(docked_mol, smiles)
            n_atoms_final = fixed_mol.GetNumAtoms()

            fixed_mol.SetProp("_Name", str(ligand_name))
            fixed_mol.SetProp("cluster_id", str(cluster_id))
            fixed_mol.SetProp("ligand_id", str(ligand_id))
            fixed_mol.SetProp("source_pdb_id", str(pdb_id))
            fixed_mol.SetProp("source_pose_rank", str(pose_rank))
            fixed_mol.SetProp("source_smiles", str(smiles))
            if pec50 is not None:
                fixed_mol.SetProp("pEC50", str(pec50))

            writer = Chem.SDWriter(str(out_ligand))
            writer.write(fixed_mol)
            writer.close()
        except Exception as e:
            lig_status = "TEMPLATE_MATCH_FAILED"
            error += f"{type(e).__name__}: {e}; "

    return {
        "cluster_id": cluster_id,
        "ligand_id": ligand_id,
        "ligand_name": ligand_name,
        "source": "test" if ligand_id.startswith("T") else "train_anchor",
        "selected_pdb_id": pdb_id,
        "selected_pose_rank": pose_rank,
        "pec50": pec50,
        "receptor_status": rec_status,
        "ligand_status": lig_status,
        "n_heavy_template": n_heavy_template,
        "n_heavy_docked": n_heavy_docked,
        "n_atoms_final": n_atoms_final,
        "error": error,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--selection", required=True,
                    help="rbfe_template_selection.csv")
    ap.add_argument("--test-full", required=True,
                    help="test_full_with_clusters_and_anchors.csv "
                         "(columns: ligand_id, Molecule Name, SMILES, "
                         "cluster_id, best_train_idx, anchor_ligand_id, "
                         "best_train_pEC50)")
    ap.add_argument("--train", required=True,
                    help="pxr-challenge_TRAIN.csv")
    ap.add_argument("--results-dir", default="results")
    ap.add_argument("--receptor-dir", default="openfe/receptors")
    ap.add_argument("--outdir", default="openfe/rbfe_inputs")
    args = ap.parse_args()

    results_dir = Path(args.results_dir)
    receptor_dir = Path(args.receptor_dir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    sel = pd.read_csv(args.selection).set_index("cluster_id")
    test_full = pd.read_csv(args.test_full)
    train = pd.read_csv(args.train).reset_index(drop=True)
    train["ligand_id"] = train.index.map(lambda i: f"A{i:04d}")

    print(f"Clusters: {len(sel)}")
    print(f"Test compounds: {len(test_full)}")
    print(f"Unique anchor A#### IDs: {test_full['anchor_ligand_id'].nunique()}")

    summary_rows = []
    n_clusters = len(sel)

    for ci, (cluster_id, sel_row) in enumerate(sel.iterrows(), 1):
        pdb_id = sel_row["selected_pdb_id"]
        pose_rank = int(sel_row["selected_pose_rank"])
        cluster_dir = outdir / str(cluster_id)
        cluster_dir.mkdir(parents=True, exist_ok=True)

        if ci % 20 == 0 or ci == 1:
            print(f"  Cluster {cluster_id} ({ci}/{n_clusters})...")

        cluster_test = test_full[test_full["cluster_id"] == cluster_id]

        # --- Test compounds ---
        for _, row in cluster_test.iterrows():
            result = extract_one(
                ligand_id=row["ligand_id"],
                ligand_name=row["Molecule Name"],
                smiles=row["SMILES"],
                pdb_id=pdb_id,
                pose_rank=pose_rank,
                cluster_id=cluster_id,
                results_dir=results_dir,
                receptor_dir=receptor_dir,
                cluster_dir=cluster_dir,
                pec50=None,
            )
            summary_rows.append(result)

        # --- All unique anchors for this cluster ---
        unique_anchor_ids = cluster_test["anchor_ligand_id"].unique()
        for anchor_lid in unique_anchor_ids:
            anchor_idx = int(anchor_lid[1:])  # A0063 -> 63
            anchor_row = train.iloc[anchor_idx]
            result = extract_one(
                ligand_id=anchor_lid,
                ligand_name=anchor_row["OCNT_ID"],
                smiles=anchor_row["SMILES"],
                pdb_id=pdb_id,
                pose_rank=pose_rank,
                cluster_id=cluster_id,
                results_dir=results_dir,
                receptor_dir=receptor_dir,
                cluster_dir=cluster_dir,
                pec50=float(anchor_row["pEC50"]),
            )
            summary_rows.append(result)

    summary = pd.DataFrame(summary_rows)
    summary_path = outdir / "extraction_summary_extended.csv"
    summary.to_csv(summary_path, index=False)

    n_ok = ((summary["receptor_status"] == "ok")
            & (summary["ligand_status"] == "ok")).sum()
    print(f"\n{n_ok} / {len(summary)} compounds fully extracted")

    failed = summary[(summary["receptor_status"] != "ok")
                      | (summary["ligand_status"] != "ok")]
    if len(failed):
        print(f"\n{len(failed)} compound(s) with issues:")
        print(failed[["cluster_id", "ligand_id", "receptor_status",
                       "ligand_status", "error"]].to_string(index=False))
    else:
        print("All extractions successful.")

    mismatch = summary[
        (summary["n_heavy_template"].notna())
        & (summary["n_heavy_docked"].notna())
        & (summary["n_heavy_template"] != summary["n_heavy_docked"])
    ]
    if len(mismatch):
        print(f"\n{len(mismatch)} heavy-atom count mismatch(es):")
        print(mismatch[["cluster_id", "ligand_id",
                         "n_heavy_template",
                         "n_heavy_docked"]].to_string(index=False))

    per_cluster = summary.groupby("cluster_id").agg(
        n_test=("ligand_id", lambda x: (x.str.startswith("T")).sum()),
        n_anchors=("ligand_id", lambda x: (x.str.startswith("A")).sum()),
        n_total=("ligand_id", "count"),
        n_ok=("ligand_status", lambda x: (x == "ok").sum()),
    ).reset_index()
    per_cluster.to_csv(outdir / "per_cluster_extraction_counts.csv",
                        index=False)
    print(f"\nCompounds per cluster: "
          f"min={per_cluster['n_total'].min()}, "
          f"max={per_cluster['n_total'].max()}, "
          f"mean={per_cluster['n_total'].mean():.1f}")
    print(f"Wrote {summary_path}")


if __name__ == "__main__":
    main()
