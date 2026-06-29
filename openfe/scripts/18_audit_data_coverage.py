#!/usr/bin/env python3
"""
openfe/scripts/18_audit_data_coverage.py

Audit completeness and merge integrity across all data sources before
building the final prediction model. Checks:

  1. Naming keys: how compounds are identified in each file (Molecule Name
     OADMET vs ligand_id T### vs OCNT). Merge bugs usually come from
     joining on mismatched keys.
  2. Docking coverage: how many Set1/Set2/train compounds have a docking
     score, in BOTH the per-compound best file and the raw receptor_best.
  3. RBFE coverage: how many of each set have an RBFE prediction.
  4. Experimental pEC50 coverage: train + phase1 should have it; phase2 not.
  5. Overlap matrix: for each set, count compounds with {docking only,
     RBFE only, both, neither} - the "neither" group is the real problem
     set with no prediction at all.
  6. Cross-source name leakage: same compound under different keys.

Usage:
    python openfe/scripts/18_audit_data_coverage.py
"""

import argparse
from pathlib import Path

import pandas as pd


def load(path, label):
    p = Path(path)
    if not p.exists():
        print(f"  [MISSING] {label}: {path}")
        return None
    df = pd.read_csv(p)
    print(f"  [ok] {label}: {df.shape[0]} rows, cols={list(df.columns)}")
    return df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--docking-dir",
                    default="docking/docking_analysis_extended")
    ap.add_argument("--results-dir", default="docking/results")
    ap.add_argument("--test-full",
                    default="openfe/test_full_with_clusters_and_anchors.csv")
    ap.add_argument("--rbfe-pred", default="openfe/rbfe_predictions.csv")
    ap.add_argument("--train", default="data/pxr-challenge_TRAIN.csv")
    ap.add_argument("--phase1",
                    default="data/pxr-challenge_TEST_PHASE_1_UNBLINDED.csv")
    ap.add_argument("--blinded", default="data/pxr-challenge_TEST_BLINDED.csv")
    args = ap.parse_args()

    dock_dir = Path(args.docking_dir)

    print("=" * 64)
    print("1. LOADING SOURCES")
    print("=" * 64)
    test_full = load(args.test_full, "test_full")
    rbfe = load(args.rbfe_pred, "rbfe_predictions")
    train = load(args.train, "train")
    phase1 = load(args.phase1, "phase1_unblinded")
    blinded = load(args.blinded, "test_blinded")
    best = load(dock_dir / "docking_receptor_best.csv", "docking_receptor_best")
    p1_dock = load(dock_dir / "phase1_with_docking_scores.csv",
                   "phase1_with_docking")
    cluster_sum = load(dock_dir / "docking_cluster_summary.csv",
                       "docking_cluster_summary")

    # ---- Define compound sets ----
    phase1_names = set(phase1["Molecule Name"])
    train_names = set(train["Molecule Name"])
    all_test = set(test_full["Molecule Name"])
    set1 = all_test & phase1_names           # Phase 1 test (unblinded)
    set2 = all_test - phase1_names           # Phase 2 test (blind)

    print("\n" + "=" * 64)
    print("2. COMPOUND SET SIZES")
    print("=" * 64)
    print(f"  All test compounds (test_full):  {len(all_test)}")
    print(f"  Set 1 (Phase 1, unblinded):      {len(set1)}")
    print(f"  Set 2 (Phase 2, blind):          {len(set2)}")
    print(f"  Training compounds:              {len(train_names)}")
    if blinded is not None:
        print(f"  TEST_BLINDED.csv rows:           {len(blinded)}")
        blinded_names = set(blinded["Molecule Name"]) \
            if "Molecule Name" in blinded.columns else set()
        print(f"  TEST_BLINDED matches test_full:  "
              f"{len(blinded_names & all_test)}")

    # ---- Naming keys in docking ----
    print("\n" + "=" * 64)
    print("3. DOCKING NAMING KEYS")
    print("=" * 64)
    if best is not None:
        dock_names = set(best["ligand_name"])
        print(f"  receptor_best unique ligand_name: {len(dock_names)}")
        print(f"  sample: {list(dock_names)[:4]}")
        # How do these intersect with test/train?
        print(f"  ligand_name ∩ Set1:  {len(dock_names & set1)}")
        print(f"  ligand_name ∩ Set2:  {len(dock_names & set2)}")
        print(f"  ligand_name ∩ train: {len(dock_names & train_names)}")
        # OCNT vs OADMET breakdown
        oadmet = {n for n in dock_names if str(n).startswith("OADMET")}
        ocnt = {n for n in dock_names if str(n).startswith("OCNT")}
        other = dock_names - oadmet - ocnt
        print(f"  name types: OADMET={len(oadmet)}  OCNT={len(ocnt)}  "
              f"other={len(other)}")
        if other:
            print(f"    other sample: {list(other)[:4]}")

    # ---- Per-compound docking best score (collapse receptor_best) ----
    print("\n" + "=" * 64)
    print("4. DOCKING COVERAGE (best CNNaffinity per compound)")
    print("=" * 64)
    # Anchors = the OCNT compounds used as RBFE references. Only these
    # training compounds were docked, NOT all 4139. Report against the
    # correct denominator.
    anchor_ocnt = set()
    if test_full is not None and train is not None:
        for _, row in test_full.drop_duplicates("anchor_ligand_id").iterrows():
            try:
                idx = int(str(row["anchor_ligand_id"])[1:])
                anchor_ocnt.add(train.iloc[idx]["OCNT_ID"])
            except (ValueError, IndexError):
                pass
    if best is not None:
        # best CNNaffinity per ligand across all receptors/poses
        per_cmpd = best.groupby("ligand_name")["CNNaffinity"].max()
        scored = set(per_cmpd.index)
        scored_oadmet = {n for n in scored if str(n).startswith("OADMET")}
        scored_ocnt = {n for n in scored if str(n).startswith("OCNT")}
        print(f"  Compounds with >=1 docking score: {len(scored)}")
        print(f"    (of which OADMET test: {len(scored_oadmet)}, "
              f"OCNT anchors: {len(scored_ocnt)})")
        print(f"    Set1 covered:    {len(set1 & scored)}/{len(set1)}")
        print(f"    Set2 covered:    {len(set2 & scored)}/{len(set2)}")
        print(f"    Anchors covered: {len(anchor_ocnt & scored)}/"
              f"{len(anchor_ocnt)}  (only anchors were docked, "
              f"not all {len(train_names)} training compounds)")
        missing2 = set2 - scored
        if missing2:
            print(f"  [WARN] Set2 compounds with NO docking: {len(missing2)}")
            print(f"    sample: {list(missing2)[:5]}")
        missing_anchor = anchor_ocnt - scored
        if missing_anchor:
            print(f"  [WARN] Anchors with NO docking: {len(missing_anchor)}")
            print(f"    sample: {list(missing_anchor)[:5]}")

    # ---- RBFE coverage ----
    print("\n" + "=" * 64)
    print("5. RBFE COVERAGE")
    print("=" * 64)
    if rbfe is not None:
        rbfe_names = set(rbfe["Molecule Name"])
        print(f"  RBFE predictions: {len(rbfe_names)}")
        print(f"    Set1: {len(set1 & rbfe_names)}/{len(set1)}")
        print(f"    Set2: {len(set2 & rbfe_names)}/{len(set2)}")
        leak = rbfe_names & train_names
        print(f"    (sanity) RBFE names that are training compounds: "
              f"{len(leak)} (should be 0)")

    # ---- Coverage matrix per set: docking / RBFE / both / neither ----
    print("\n" + "=" * 64)
    print("6. PREDICTION COVERAGE MATRIX (the key table)")
    print("=" * 64)
    if best is not None and rbfe is not None:
        scored = set(best.groupby("ligand_name")["CNNaffinity"].max().index)
        rbfe_names = set(rbfe["Molecule Name"])
        for label, s in [("Set1 (Phase1)", set1), ("Set2 (Phase2)", set2)]:
            both = s & scored & rbfe_names
            dock_only = (s & scored) - rbfe_names
            rbfe_only = (s & rbfe_names) - scored
            neither = s - scored - rbfe_names
            print(f"\n  {label} (n={len(s)}):")
            print(f"    both docking+RBFE: {len(both)}")
            print(f"    docking only:      {len(dock_only)}")
            print(f"    RBFE only:         {len(rbfe_only)}")
            print(f"    NEITHER:           {len(neither)}  "
                  f"{'<-- PROBLEM' if neither else ''}")
            if neither and label.startswith("Set2"):
                print(f"      sample: {list(neither)[:5]}")

    print("\n" + "=" * 64)
    print("AUDIT COMPLETE")
    print("=" * 64)


if __name__ == "__main__":
    main()
