#!/usr/bin/env python3
"""
analysis/03_cluster_to_crystal_template.py

For each test-set cluster (from 01_cluster_test_set.py), find the
crystal structure (from the 64-structure re-refined PXR ensemble, via
analyze_pxr_ensemble.py's pxr_structure_inventory.csv) whose bound
ligand is most chemically similar to any compound in that cluster.

This was used to check whether any of the 64 PXR crystal ligands could
serve as a pose/template donor for the test-set chemotype. RESULT (as
of the original run): best similarities were ~0.12-0.21 across all 124
clusters - i.e. NO good template exists in the crystal ensemble. This
script is kept for re-use (e.g. if the 184-structure release adds new
ligand chemotypes, or for sanity-checking future additions).

Produces:
  - cluster_template_mapping.csv : one row per cluster with the best-
    matching crystal structure, its ligand identity, and the similarity
    score.

Usage:
    pip install rdkit pandas numpy
    python 03_cluster_to_crystal_template.py \
        --test-clusters test_with_clusters.csv \
        --structure-inventory pxr_structure_inventory.csv
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")


def get_fp(smi, radius=2, nbits=2048):
    if pd.isna(smi):
        return None
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    return AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=nbits)


def clean_pdb_smiles(smi):
    """PDB-connectivity-derived SMILES often have spurious explicit Hs
    and ambiguous bond orders. Strip Hs and re-sanitize loosely so we
    can at least compute a Morgan fingerprint for rough similarity
    matching (NOT a validated structure)."""
    if pd.isna(smi):
        return None
    mol = Chem.MolFromSmiles(smi, sanitize=False)
    if mol is None:
        return None
    try:
        mol = Chem.RemoveHs(mol, sanitize=False)
        Chem.SanitizeMol(
            mol,
            sanitizeOps=Chem.SANITIZE_ALL
            ^ Chem.SANITIZE_PROPERTIES
            ^ Chem.SANITIZE_ADJUSTHS,
        )
        return Chem.MolToSmiles(mol)
    except Exception:
        return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--test-clusters", required=True,
                    help="Output of 01_cluster_test_set.py "
                         "(test_with_clusters.csv)")
    ap.add_argument("--structure-inventory", required=True,
                    help="Output of analyze_pxr_ensemble.py "
                         "(pxr_structure_inventory.csv)")
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    test = pd.read_csv(args.test_clusters)
    inv = pd.read_csv(args.structure_inventory)
    print(f"Test compounds: {len(test)}, clusters: {test['cluster_id'].nunique()}")
    print(f"Crystal structures: {len(inv)}")

    inv["clean_smiles"] = inv["ligand_smiles"].apply(clean_pdb_smiles)
    inv["fp"] = inv["clean_smiles"].apply(get_fp)
    xtal_valid = inv[inv["fp"].notna()].reset_index(drop=True)
    print(f"Crystal ligands with usable fingerprints: "
          f"{len(xtal_valid)} / {len(inv)}")

    test["fp"] = test["SMILES"].apply(get_fp)

    xtal_fps = xtal_valid["fp"].tolist()
    xtal_ids = xtal_valid["pdb_id"].tolist()
    xtal_ligs = xtal_valid["ligand_resname"].tolist()

    rows = {}
    for cid, group in test.groupby("cluster_id"):
        best_sim, best_pdb, best_lig, best_idx = -1, None, None, None
        for idx, row in group.iterrows():
            sims = np.array(DataStructs.BulkTanimotoSimilarity(row["fp"], xtal_fps))
            j = int(sims.argmax())
            if sims[j] > best_sim:
                best_sim, best_pdb, best_lig, best_idx = sims[j], xtal_ids[j], xtal_ligs[j], idx
        rows[cid] = {
            "cluster_size": len(group),
            "best_sim_to_xtal": round(best_sim, 3),
            "best_pdb_template": best_pdb,
            "template_ligand": best_lig,
            "representative_test_compound": test.loc[best_idx, "Molecule Name"],
            "representative_smiles": test.loc[best_idx, "SMILES"],
        }

    ct = pd.DataFrame.from_dict(rows, orient="index").sort_values(
        "cluster_size", ascending=False)

    print("\nBest-template similarity distribution across clusters:")
    print(ct["best_sim_to_xtal"].describe())
    for thr in (0.2, 0.3, 0.4):
        n_clusters = (ct["best_sim_to_xtal"] >= thr).sum()
        n_compounds = ct[ct["best_sim_to_xtal"] >= thr]["cluster_size"].sum()
        print(f"  clusters with best_sim >= {thr}: {n_clusters} "
              f"({n_compounds} compounds)")

    out_path = outdir / "cluster_template_mapping.csv"
    ct.to_csv(out_path)
    print(f"\nWrote {out_path}")


if __name__ == "__main__":
    main()
