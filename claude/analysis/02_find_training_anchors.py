#!/usr/bin/env python3
"""
analysis/02_find_training_anchors.py

For each compound in the 513-compound test set, find its nearest
neighbor (by ECFP4 Tanimoto similarity) in the 4,139-compound training
set. This identifies candidate calibration/RBFE-network anchors with
known pEC50.

Produces:
  - test_with_train_anchors.csv : test set + columns
        max_train_sim       (Tanimoto similarity to nearest train analog)
        best_train_idx      (row index into the training set)
        best_train_pEC50    (pEC50 of that nearest train analog)
        best_train_OCNT     (OCNT_ID of that nearest train analog, if
                              present in the training data)

Usage:
    pip install rdkit pandas numpy
    python 02_find_training_anchors.py \
        --test pxr-challenge_TEST_BLINDED.csv \
        --train pxr-challenge_TRAIN.csv
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--test", required=True)
    ap.add_argument("--train", required=True)
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    test = pd.read_csv(args.test)
    train = pd.read_csv(args.train)
    print(f"Test: {len(test)}, Train: {len(train)}")

    test["fp"] = test["SMILES"].apply(get_fp)
    train["fp"] = train["SMILES"].apply(get_fp)

    valid_test = test[test["fp"].notna()].reset_index(drop=True)
    valid_train = train[train["fp"].notna()].reset_index(drop=True)
    print(f"Valid fps - test: {len(valid_test)}, train: {len(valid_train)}")

    train_fps = valid_train["fp"].tolist()
    test_fps = valid_test["fp"].tolist()

    max_sims, best_idx = [], []
    for fp in test_fps:
        sims = np.array(DataStructs.BulkTanimotoSimilarity(fp, train_fps))
        j = int(sims.argmax())
        max_sims.append(float(sims[j]))
        best_idx.append(j)

    valid_test["max_train_sim"] = max_sims
    valid_test["best_train_idx"] = best_idx
    valid_test["best_train_pEC50"] = valid_train.iloc[best_idx]["pEC50"].values

    # OCNT_ID may or may not be present depending on training CSV version
    if "OCNT_ID" in valid_train.columns:
        valid_test["best_train_OCNT"] = valid_train.iloc[best_idx]["OCNT_ID"].values

    print("\nDistribution of max similarity (test -> nearest train):")
    print(valid_test["max_train_sim"].describe())

    for thr in (0.4, 0.5, 0.6, 0.7):
        n = (valid_test["max_train_sim"] >= thr).sum()
        print(f"  test compounds with a train analog >= {thr}: {n}")

    out_path = outdir / "test_with_train_anchors.csv"
    valid_test.drop(columns=["fp"]).to_csv(out_path, index=False)
    print(f"\nWrote {out_path}")


if __name__ == "__main__":
    main()
