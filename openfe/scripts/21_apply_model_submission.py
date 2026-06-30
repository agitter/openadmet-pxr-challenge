#!/usr/bin/env python3
"""
openfe/scripts/21_apply_model_submission.py

Apply the frozen model (chosen_model.json) to ALL 513 test compounds and
write the competition submission file.

Pure application - NO fitting. Loads the frozen calibrations + blend
weight, applies:
    pred = w*rbfe_cal + (1-w)*dock_cal   where RBFE present
    pred = dock_cal                       where RBFE absent
with dock_cal = dock_slope*CNNaffinity + dock_intercept and
     rbfe_cal = rbfe_slope*pred_pEC50_raw + rbfe_intercept.

Every test compound has a docking score (verified in audit), so there
are no missing predictions - the evaluator requires no NaNs.

Output: submission.csv with columns [SMILES, Molecule Name, pEC50].

Usage:
    python openfe/scripts/21_apply_model_submission.py \
        --model openfe/chosen_model.json \
        --blinded data/pxr-challenge_TEST_BLINDED.csv \
        --rbfe-pred openfe/rbfe_predictions.csv \
        --receptor-best docking/docking_analysis_extended/docking_receptor_best.csv \
        --phase1 data/pxr-challenge_TEST_PHASE_1_UNBLINDED.csv \
        --outdir openfe
"""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", default="openfe/chosen_model.json")
    ap.add_argument("--blinded", default="data/pxr-challenge_TEST_BLINDED.csv")
    ap.add_argument("--rbfe-pred", default="openfe/rbfe_predictions.csv")
    ap.add_argument("--receptor-best",
                    default="docking/docking_analysis_extended/"
                            "docking_receptor_best.csv")
    ap.add_argument("--phase1",
                    default="data/pxr-challenge_TEST_PHASE_1_UNBLINDED.csv")
    ap.add_argument("--test-full",
                    default="openfe/test_full_with_clusters_and_anchors.csv")
    ap.add_argument("--outdir", default="openfe")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    with open(args.model) as f:
        model = json.load(f)
    cal = model["calibrations"]
    w = model["winner_params"]["w"]
    print(f"Loaded model: {model['winner']} (w={w})")
    print(f"  dock_cal: {cal['dock_slope']:.4f}*CNNaffinity "
          f"+ {cal['dock_intercept']:.4f}")
    print(f"  rbfe_cal: {cal['rbfe_slope']:.4f}*pred_pEC50_raw "
          f"+ {cal['rbfe_intercept']:.4f}")

    blinded = pd.read_csv(args.blinded)        # all 513, the submission set
    rbfe = pd.read_csv(args.rbfe_pred)
    best = pd.read_csv(args.receptor_best)
    phase1 = pd.read_csv(args.phase1)

    # Best CNNaffinity per compound (all 513 covered)
    dock = (best.sort_values("CNNaffinity", ascending=False)
                .groupby("ligand_name", as_index=False).first()
                [["ligand_name", "CNNaffinity"]]
                .rename(columns={"ligand_name": "Molecule Name"}))

    # Assemble per-compound inputs for all 513 (carry SMILES from blinded)
    keep = ["Molecule Name"]
    if "SMILES" in blinded.columns:
        keep.append("SMILES")
    sub = blinded[keep].copy()
    sub = sub.merge(dock, on="Molecule Name", how="left")
    sub = sub.merge(rbfe[["Molecule Name", "pred_pEC50_raw"]],
                    on="Molecule Name", how="left")
    sub["has_rbfe"] = sub["pred_pEC50_raw"].notna()

    # Sanity: every compound must have docking
    missing_dock = sub["CNNaffinity"].isna().sum()
    if missing_dock:
        raise SystemExit(f"ERROR: {missing_dock} compounds lack docking "
                         f"scores; cannot produce a complete submission.")

    # Apply frozen model
    dock_cal = cal["dock_slope"] * sub["CNNaffinity"].values \
        + cal["dock_intercept"]
    rbfe_cal = np.where(
        sub["has_rbfe"].values,
        cal["rbfe_slope"] * sub["pred_pEC50_raw"].values + cal["rbfe_intercept"],
        cal["train_mean"])
    wv = np.where(sub["has_rbfe"].values, w, 0.0)
    sub["pEC50"] = wv * rbfe_cal + (1 - wv) * dock_cal

    # Tag set membership + prediction source for our own records
    phase1_names = set(phase1["Molecule Name"])
    sub["set"] = sub["Molecule Name"].apply(
        lambda x: "set1_phase1" if x in phase1_names else "set2_phase2")
    sub["source"] = np.where(sub["has_rbfe"], "blend_w0.5", "docking_only")

    # ---- Submission file: SMILES + Molecule Name + pEC50, all 513 ----
    # (validator requires SMILES, Molecule Name as IDs; pEC50 as value;
    #  no nulls, no duplicates, all finite, exactly 513 rows)
    out_cols = ["Molecule Name", "pEC50"]
    if "SMILES" in sub.columns:
        out_cols = ["SMILES", "Molecule Name", "pEC50"]
    submission = sub[out_cols].copy()

    # Validation mirroring the organizers' activity_validation.py
    problems = []
    if submission["pEC50"].isna().any():
        problems.append("NaN in pEC50")
    if not np.isfinite(submission["pEC50"].to_numpy()).all():
        problems.append("non-finite pEC50")
    if submission["Molecule Name"].duplicated().any():
        problems.append("duplicate Molecule Name")
    if "SMILES" in submission.columns and submission["SMILES"].isna().any():
        problems.append("missing SMILES")
    if submission["Molecule Name"].isna().any():
        problems.append("missing Molecule Name")
    if len(submission) != 513:
        problems.append(f"row count {len(submission)} != 513")
    if problems:
        raise SystemExit("Submission validation FAILED: " + "; ".join(problems))
    print("Submission validation: PASSED "
          "(513 rows, no NaN/inf, no duplicates, SMILES present)")

    sub_path = outdir / "submission.csv"
    submission.to_csv(sub_path, index=False)

    # ---- Detailed version for our records ----
    detail_path = outdir / "submission_detailed.csv"
    sub[["Molecule Name", "set", "source", "CNNaffinity",
         "pred_pEC50_raw", "has_rbfe", "pEC50"]].to_csv(detail_path,
                                                        index=False)

    # ---- Visualizations of the submitted predictions ----
    # Bring in Phase 1 truth for the calibration-quality panel
    truth = dict(zip(phase1["Molecule Name"], phase1["pEC50"]))
    sub["exp_pEC50"] = sub["Molecule Name"].map(truth)

    fig, axes = plt.subplots(2, 3, figsize=(17, 10))
    set_colors = {"set1_phase1": "#fd8d3c", "set2_phase2": "#74c476"}

    # Panel 1: submitted pEC50 distribution, Set1 vs Set2
    ax = axes[0, 0]
    for s, c in set_colors.items():
        vals = sub[sub["set"] == s]["pEC50"]
        ax.hist(vals, bins=30, alpha=0.55, label=f"{s} (n={len(vals)})",
                color=c, density=True)
    ax.set_xlabel("Submitted pEC50"); ax.set_ylabel("Density")
    ax.set_title("Submitted pEC50 by set"); ax.legend(fontsize=8)

    # Panel 2: submitted pEC50 by prediction source
    ax = axes[0, 1]
    for src, c in [("blend_w0.5", "#3182bd"), ("docking_only", "#969696")]:
        vals = sub[sub["source"] == src]["pEC50"]
        ax.hist(vals, bins=30, alpha=0.55, label=f"{src} (n={len(vals)})",
                color=c, density=True)
    ax.set_xlabel("Submitted pEC50"); ax.set_ylabel("Density")
    ax.set_title("Submitted pEC50 by source"); ax.legend(fontsize=8)

    # Panel 3: submitted vs experimental (Set 1, where truth known)
    ax = axes[0, 2]
    s1 = sub[sub["set"] == "set1_phase1"].dropna(subset=["exp_pEC50"])
    cov = s1["has_rbfe"].values
    ax.scatter(s1["pEC50"].values[~cov], s1["exp_pEC50"].values[~cov],
               s=16, alpha=0.5, c="#969696", label="docking-only")
    ax.scatter(s1["pEC50"].values[cov], s1["exp_pEC50"].values[cov],
               s=16, alpha=0.5, c="#3182bd", label="blend")
    lims = [min(s1["pEC50"].min(), s1["exp_pEC50"].min()) - 0.3,
            max(s1["pEC50"].max(), s1["exp_pEC50"].max()) + 0.3]
    ax.plot(lims, lims, "k--", lw=0.7, alpha=0.6)
    ax.set_xlim(lims); ax.set_ylim(lims)
    ax.set_xlabel("Submitted pEC50"); ax.set_ylabel("Experimental pEC50")
    ax.set_title("Set 1 check: submitted vs experimental")
    ax.legend(fontsize=8)

    # Panel 4: predicted pEC50 vs CNNaffinity (shows docking's role)
    ax = axes[1, 0]
    for s, c in set_colors.items():
        d = sub[sub["set"] == s]
        ax.scatter(d["CNNaffinity"], d["pEC50"], s=12, alpha=0.4,
                   c=c, label=s)
    ax.set_xlabel("Docking CNNaffinity"); ax.set_ylabel("Submitted pEC50")
    ax.set_title("Submitted pEC50 vs CNNaffinity"); ax.legend(fontsize=8)

    # Panel 5: raw RBFE vs submitted (only blend compounds) - shows how much
    # the blend pulled RBFE toward the calibrated scale
    ax = axes[1, 1]
    bl = sub[sub["has_rbfe"]]
    ax.scatter(bl["pred_pEC50_raw"], bl["pEC50"], s=14, alpha=0.5,
               c="#3182bd")
    ax.set_xlabel("Raw RBFE pEC50 (uncalibrated)")
    ax.set_ylabel("Submitted pEC50 (blended)")
    ax.set_title("Blend compounds: raw RBFE -> submitted")

    # Panel 6: ECDF of submitted Set 2 vs known training/Phase1 pEC50
    ax = axes[1, 2]
    def ecdf(v):
        v = np.sort(v); return v, np.arange(1, len(v) + 1) / len(v)
    x2, y2 = ecdf(sub[sub["set"] == "set2_phase2"]["pEC50"].values)
    ax.plot(x2, y2, label="Set 2 submitted", color="#74c476", lw=2)
    xp, yp = ecdf(phase1["pEC50"].values)
    ax.plot(xp, yp, label="Phase 1 experimental", color="#fd8d3c",
            lw=2, ls="--")
    ax.set_xlabel("pEC50"); ax.set_ylabel("Cumulative fraction")
    ax.set_title("Set 2 submitted vs Phase 1 experimental (ECDF)")
    ax.legend(fontsize=8)

    plt.tight_layout()
    viz_path = outdir / "submission_visualizations.png"
    plt.savefig(viz_path, dpi=130, bbox_inches="tight")

    print(f"\nSubmission: {len(submission)} compounds")
    print(f"  Set 1 (Phase 1): {(sub['set']=='set1_phase1').sum()}")
    print(f"  Set 2 (Phase 2): {(sub['set']=='set2_phase2').sum()}")
    print(f"  blend (has RBFE): {(sub['source']=='blend_w0.5').sum()}")
    print(f"  docking only:     {(sub['source']=='docking_only').sum()}")
    print(f"\nPredicted pEC50 distribution:")
    print(f"  mean={submission['pEC50'].mean():.2f}  "
          f"median={submission['pEC50'].median():.2f}  "
          f"min={submission['pEC50'].min():.2f}  "
          f"max={submission['pEC50'].max():.2f}")
    # Set 2 only (the scored set)
    s2 = sub[sub["set"] == "set2_phase2"]["pEC50"]
    print(f"\nSet 2 (scored) pEC50: mean={s2.mean():.2f} "
          f"min={s2.min():.2f} max={s2.max():.2f}")
    print(f"\nWrote {sub_path}")
    print(f"Wrote {detail_path}")
    print(f"Wrote {viz_path}")


if __name__ == "__main__":
    main()
