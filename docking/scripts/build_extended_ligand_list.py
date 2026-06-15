#!/usr/bin/env python3
"""
docking/scripts/build_extended_ligand_list.py

Build the full ligand list for the extended single-cluster-receptor-set
docking campaign: all 513 test compounds + all unique nearest-train-
analog "anchors" (from 02_find_training_anchors.py's output), each to
be docked against ALL 62 receptors (same as the original 124-cluster-
representative campaign).

This extended campaign exists to enable a post-hoc RBFE-vs-docking
disagreement check: for any RBFE edge (A -> B), having CNNaffinity /
minimizedAffinity for both A and B (across the same receptor set) lets
us compute a cheap docking-based delta-delta-G estimate to compare
against the RBFE-computed value.

Ligand IDs (used as work-unit / results directory names, analogous to
cluster_id in the original campaign):
  Test compounds:   T<3-digit zero-padded row index into TEST_BLINDED>
                     e.g. T000 .. T512
  Training anchors: A<4-digit zero-padded row index into TRAIN>
                     e.g. A0063 .. A3223

These prefixes ensure no collision with the original campaign's purely
numeric cluster_id directories (e.g. "21", "401") already present under
results/.

Inputs:
  --test-with-anchors  test_with_train_anchors.csv (output of
                        02_find_training_anchors.py): must have
                        'Molecule Name', 'SMILES', 'best_train_idx'
  --train              pxr-challenge_TRAIN.csv: must have 'OCNT_ID',
                        'SMILES'

Output:
  all_rbfe_ligands.csv: ligand_id, ligand_name, smiles, source
    source in {test, train_anchor}

Usage:
    pip install pandas
    python build_extended_ligand_list.py \
        --test-with-anchors test_with_train_anchors.csv \
        --train pxr-challenge_TRAIN.csv \
        --outdir docking
"""

import argparse
from pathlib import Path

import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--test-with-anchors", required=True)
    ap.add_argument("--train", required=True)
    ap.add_argument("--outdir", default="docking")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    test = pd.read_csv(args.test_with_anchors)
    train = pd.read_csv(args.train)

    print(f"Test compounds: {len(test)}")

    # --- Test compounds: T<3-digit row index> ---
    test_rows = []
    for idx, row in test.reset_index(drop=True).iterrows():
        test_rows.append({
            "ligand_id": f"T{idx:03d}",
            "ligand_name": row["Molecule Name"],
            "smiles": row["SMILES"],
            "source": "test",
        })

    # --- Unique training anchors: A<4-digit row index into TRAIN> ---
    unique_anchor_idx = sorted(test["best_train_idx"].unique())
    print(f"Unique training anchors: {len(unique_anchor_idx)}")

    anchor_rows = []
    for idx in unique_anchor_idx:
        idx = int(idx)
        row = train.iloc[idx]
        anchor_rows.append({
            "ligand_id": f"A{idx:04d}",
            "ligand_name": row["OCNT_ID"],
            "smiles": row["SMILES"],
            "source": "train_anchor",
        })

    combined = pd.DataFrame(test_rows + anchor_rows)
    out_path = outdir / "all_rbfe_ligands.csv"
    combined.to_csv(out_path, index=False)

    print(f"\nTotal ligands: {len(combined)} "
          f"({len(test_rows)} test + {len(anchor_rows)} anchors)")
    print(f"Wrote {out_path}")
    print(f"\nWith 62 receptors: {len(combined)*62} (ligand, receptor) pairs "
          f"-> {len(combined)} HTCondor jobs (one ligand x all receptors "
          f"each, same pattern as the original 124-job campaign)")


if __name__ == "__main__":
    main()
