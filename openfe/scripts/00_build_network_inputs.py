#!/usr/bin/env python3
"""
openfe/scripts/00_build_network_inputs.py

Build the per-compound mapping table needed for RBFE network design:
  - Joins the 513 test compounds with their cluster assignments
    (from 01_cluster_test_set.py output) and nearest training anchors
    (from 02_find_training_anchors.py output)
  - Assigns T### ligand IDs (matching the extended docking campaign)
  - Assigns A#### anchor ligand IDs (matching the extended docking campaign)

This table is the primary input to 03_extract_all_rbfe_inputs.py.

Inputs:
  --test        pxr-challenge_TEST_BLINDED.csv
  --clusters    test_with_clusters.csv (output of analysis/01_cluster_test_set.py)
  --anchors     test_with_train_anchors.csv (output of analysis/02_find_training_anchors.py)

Output:
  --outdir/test_full_with_clusters_and_anchors.csv
      columns: ligand_id, Molecule Name, SMILES, cluster_id,
               best_train_idx, anchor_ligand_id, best_train_pEC50

Usage:
    python openfe/scripts/00_build_network_inputs.py \
        --test data/pxr-challenge_TEST_BLINDED.csv \
        --clusters analysis/outputs/test_with_clusters.csv \
        --anchors analysis/outputs/test_with_train_anchors.csv \
        --outdir openfe
"""

import argparse
from pathlib import Path

import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--test", required=True,
                    help="pxr-challenge_TEST_BLINDED.csv")
    ap.add_argument("--clusters", required=True,
                    help="test_with_clusters.csv "
                         "(output of analysis/01_cluster_test_set.py)")
    ap.add_argument("--anchors", required=True,
                    help="test_with_train_anchors.csv "
                         "(output of analysis/02_find_training_anchors.py)")
    ap.add_argument("--outdir", default="openfe",
                    help="Output directory (default: openfe/)")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    test = pd.read_csv(args.test).reset_index(drop=True)
    test["ligand_id"] = test.index.map(lambda i: f"T{i:03d}")

    clusters = pd.read_csv(args.clusters)
    anchors = pd.read_csv(args.anchors)

    print(f"Test compounds: {len(test)}")
    print(f"Cluster assignments: {len(clusters)} "
          f"({clusters['cluster_id'].nunique()} unique clusters)")
    print(f"Anchor assignments: {len(anchors)}")

    test_full = test.merge(
        clusters[["Molecule Name", "cluster_id"]], on="Molecule Name")
    test_full = test_full.merge(
        anchors[["Molecule Name", "best_train_idx", "best_train_pEC50"]],
        on="Molecule Name")
    test_full["anchor_ligand_id"] = test_full["best_train_idx"].map(
        lambda i: f"A{i:04d}")

    n_unique_anchors = test_full["anchor_ligand_id"].nunique()
    print(f"\nJoined table: {len(test_full)} rows")
    print(f"Unique anchor A#### IDs: {n_unique_anchors}")
    print(f"Clusters with >1 unique anchor: "
          f"{(test_full.groupby('cluster_id')['anchor_ligand_id'].nunique() > 1).sum()}")

    out_path = outdir / "test_full_with_clusters_and_anchors.csv"
    test_full.to_csv(out_path, index=False)
    print(f"\nWrote {out_path}")


if __name__ == "__main__":
    main()
