#!/usr/bin/env python3
"""
openfe/scripts/17_compare_methods.py

Head-to-head comparison of prediction methods against experimental pEC50,
with diagnostic visualizations. Compares:
  - Docking scores (CNNaffinity, CNNscore, minimizedAffinity)
  - RBFE-propagated pEC50
against known experimental pEC50, across compound categories
(training, Phase 1 test, Phase 2 test).

Goal: decide which signal (or combination) to trust for the Phase 2
submission, and document it for the writeup.

Produces:
  - openfe/method_comparison_metrics.csv   (correlations per method/subset)
  - openfe/method_comparison.png           (multi-panel figure)
  - console summary

Usage:
    python openfe/scripts/17_compare_methods.py \
        --phase1-docking docking/docking_analysis_extended/phase1_with_docking_scores.csv \
        --train-docking docking/docking_analysis_extended/train_anchors_with_docking_scores.csv \
        --receptor-best docking/docking_analysis_extended/docking_receptor_best.csv \
        --rbfe-pred openfe/rbfe_predictions.csv \
        --edge-convergence openfe/edge_convergence.csv \
        --test-full openfe/test_full_with_clusters_and_anchors.csv \
        --train data/pxr-challenge_TRAIN.csv \
        --phase1 data/pxr-challenge_TEST_PHASE_1_UNBLINDED.csv \
        --outdir openfe
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr


def corr_block(x, y):
    """Return (n, pearson_r, spearman_rho, mae) for paired arrays,
    dropping NaNs."""
    mask = ~(pd.isna(x) | pd.isna(y))
    x, y = np.asarray(x)[mask], np.asarray(y)[mask]
    if len(x) < 5:
        return len(x), np.nan, np.nan, np.nan
    r, _ = pearsonr(x, y)
    rho, _ = spearmanr(x, y)
    mae = np.mean(np.abs(x - y))
    return len(x), r, rho, mae


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--phase1-docking",
                    default="docking/docking_analysis_extended/"
                            "phase1_with_docking_scores.csv")
    ap.add_argument("--train-docking",
                    default="docking/docking_analysis_extended/"
                            "train_anchors_with_docking_scores.csv")
    ap.add_argument("--receptor-best",
                    default="docking/docking_analysis_extended/"
                            "docking_receptor_best.csv")
    ap.add_argument("--rbfe-pred", default="openfe/rbfe_predictions.csv")
    ap.add_argument("--edge-convergence",
                    default="openfe/edge_convergence.csv")
    ap.add_argument("--test-full",
                    default="openfe/test_full_with_clusters_and_anchors.csv")
    ap.add_argument("--train", default="data/pxr-challenge_TRAIN.csv")
    ap.add_argument("--phase1",
                    default="data/pxr-challenge_TEST_PHASE_1_UNBLINDED.csv")
    ap.add_argument("--outdir", default="openfe")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    p1_dock = pd.read_csv(args.phase1_docking)
    rbfe = pd.read_csv(args.rbfe_pred)
    test_full = pd.read_csv(args.test_full)
    train = pd.read_csv(args.train)
    phase1 = pd.read_csv(args.phase1)

    # ---- Assemble a master per-compound table ----
    # Experimental pEC50 source: train + phase1 (both known). Phase 2 unknown.
    # Key by BOTH Molecule Name (OADMET) and OCNT_ID, because docking files
    # key anchors by OCNT_ID while test/train files use OADMET Molecule Name.
    exp = {}
    train_ocnt_names = set()
    for _, r in train.iterrows():
        exp[r["Molecule Name"]] = r["pEC50"]
        if "OCNT_ID" in train.columns and pd.notna(r.get("OCNT_ID")):
            exp[r["OCNT_ID"]] = r["pEC50"]
            train_ocnt_names.add(r["OCNT_ID"])
    for _, r in phase1.iterrows():
        exp[r["Molecule Name"]] = r["pEC50"]

    phase1_names = set(phase1["Molecule Name"])
    train_names = set(train["Molecule Name"])

    def category(name):
        if name in train_names or name in train_ocnt_names:
            return "train"
        if name in phase1_names:
            return "phase1"
        return "phase2"

    # Docking per-compound: use the FULL receptor_best file (all 513 test
    # compounds + 89 anchors), collapsing to best CNNaffinity per ligand.
    # The phase1-only file would miss Set2 and anchors.
    best = pd.read_csv(args.receptor_best)
    dock = (best.sort_values("CNNaffinity", ascending=False)
                .groupby("ligand_name", as_index=False)
                .first()[["ligand_name", "CNNaffinity", "CNNscore",
                          "minimizedAffinity"]])
    dock = dock.rename(columns={"ligand_name": "Molecule Name"})

    # RBFE predictions
    rbfe_small = rbfe[["Molecule Name", "pred_pEC50_raw", "n_hops"]].copy()

    # Master merge
    master = pd.DataFrame({"Molecule Name": list(
        set(dock["Molecule Name"]) | set(rbfe_small["Molecule Name"])
        | set(exp.keys()))})
    master["exp_pEC50"] = master["Molecule Name"].map(exp)
    master["category"] = master["Molecule Name"].apply(category)
    master = master.merge(dock, on="Molecule Name", how="left")
    master = master.merge(rbfe_small, on="Molecule Name", how="left")
    master.to_csv(outdir / "method_comparison_master.csv", index=False)

    # ---- Metrics: each predictor vs experimental pEC50, by subset ----
    # "train" subset is meaningless (only 89 anchors were docked, and RBFE
    # has no training predictions). Report anchors (docked training compounds)
    # separately so the denominator is correct.
    anchor_names = set(master[(master["category"] == "train") &
                              master["CNNaffinity"].notna()]["Molecule Name"])
    metric_rows = []
    predictors = ["CNNaffinity", "CNNscore", "minimizedAffinity",
                  "pred_pEC50_raw"]
    for subset_name, subset in [
        ("phase1", master[master["category"] == "phase1"]),
        ("anchors_docked",
         master[master["Molecule Name"].isin(anchor_names)]),
    ]:
        for pred in predictors:
            if pred not in subset.columns:
                continue
            n, r, rho, mae = corr_block(subset[pred], subset["exp_pEC50"])
            metric_rows.append({
                "subset": subset_name, "predictor": pred,
                "n": n, "pearson": r, "spearman": rho,
                "mae_vs_pEC50": mae if pred == "pred_pEC50_raw" else np.nan,
            })

    metrics = pd.DataFrame(metric_rows)
    metrics.to_csv(outdir / "method_comparison_metrics.csv", index=False)

    print("=" * 64)
    print("METHOD vs EXPERIMENTAL pEC50")
    print("=" * 64)
    print(metrics.to_string(index=False))

    # ---- Visualizations ----
    fig, axes = plt.subplots(2, 3, figsize=(16, 10))

    cat_colors = {"train": "#9ecae1", "phase1": "#fd8d3c",
                  "phase2": "#74c476"}

    # Panel 1: CNNaffinity vs exp pEC50 (phase1)
    ax = axes[0, 0]
    p1 = master[master["category"] == "phase1"]
    ax.scatter(p1["CNNaffinity"], p1["exp_pEC50"], alpha=0.5,
               c=cat_colors["phase1"], s=20)
    n, r, rho, _ = corr_block(p1["CNNaffinity"], p1["exp_pEC50"])
    ax.set_xlabel("Docking CNNaffinity")
    ax.set_ylabel("Experimental pEC50")
    ax.set_title(f"Docking CNNaffinity vs pEC50 (Phase 1)\n"
                 f"r={r:.3f}  rho={rho:.3f}  n={n}")

    # Panel 2: RBFE pred vs exp pEC50 (phase1)
    ax = axes[0, 1]
    ax.scatter(p1["pred_pEC50_raw"], p1["exp_pEC50"], alpha=0.5,
               c=cat_colors["phase1"], s=20)
    n, r, rho, mae = corr_block(p1["pred_pEC50_raw"], p1["exp_pEC50"])
    lims = [0, 10]
    ax.plot(lims, lims, "k--", lw=0.8, alpha=0.5)
    ax.set_xlabel("RBFE predicted pEC50 (raw)")
    ax.set_ylabel("Experimental pEC50")
    ax.set_title(f"RBFE vs pEC50 (Phase 1)\n"
                 f"r={r:.3f}  rho={rho:.3f}  MAE={mae:.2f}  n={n}")
    ax.set_xlim(-6, 16)

    # Panel 3: docking vs RBFE (do they agree? phase1)
    ax = axes[0, 2]
    sc = ax.scatter(p1["CNNaffinity"], p1["pred_pEC50_raw"],
                    c=p1["exp_pEC50"], cmap="viridis", alpha=0.6, s=20)
    n, r, rho, _ = corr_block(p1["CNNaffinity"], p1["pred_pEC50_raw"])
    ax.set_xlabel("Docking CNNaffinity")
    ax.set_ylabel("RBFE predicted pEC50")
    ax.set_title(f"Docking vs RBFE (Phase 1)\n"
                 f"agreement r={r:.3f}  rho={rho:.3f}")
    plt.colorbar(sc, ax=ax, label="exp pEC50")

    # Panel 4: experimental pEC50 distributions by category
    ax = axes[1, 0]
    for cat in ["train", "phase1"]:
        vals = master[master["category"] == cat]["exp_pEC50"].dropna()
        ax.hist(vals, bins=30, alpha=0.5, label=f"{cat} (n={len(vals)})",
                color=cat_colors[cat], density=True)
    ax.set_xlabel("Experimental pEC50")
    ax.set_ylabel("Density")
    ax.set_title("Experimental pEC50 distribution\n(narrow = hard to rank)")
    ax.legend()

    # Panel 5: CNNaffinity distribution by category (incl phase2)
    ax = axes[1, 1]
    for cat in ["train", "phase1", "phase2"]:
        vals = master[master["category"] == cat]["CNNaffinity"].dropna()
        if len(vals):
            ax.hist(vals, bins=30, alpha=0.5,
                    label=f"{cat} (n={len(vals)})",
                    color=cat_colors[cat], density=True)
    ax.set_xlabel("Docking CNNaffinity")
    ax.set_ylabel("Density")
    ax.set_title("CNNaffinity distribution by set")
    ax.legend()

    # Panel 6: RBFE error vs |prediction deviation| (does high error = bad?)
    ax = axes[1, 2]
    p1v = p1.dropna(subset=["pred_pEC50_raw", "exp_pEC50"]).copy()
    p1v["abs_err"] = (p1v["pred_pEC50_raw"] - p1v["exp_pEC50"]).abs()
    ax.scatter(p1v["n_hops"], p1v["abs_err"], alpha=0.5,
               c=cat_colors["phase1"], s=20)
    n, r, rho, _ = corr_block(p1v["n_hops"], p1v["abs_err"])
    ax.set_xlabel("RBFE path length (hops)")
    ax.set_ylabel("|RBFE pred - exp pEC50|")
    ax.set_title(f"Does longer path = worse?\n"
                 f"r={r:.3f}  rho={rho:.3f}")

    plt.tight_layout()
    fig_path = outdir / "method_comparison.png"
    plt.savefig(fig_path, dpi=130, bbox_inches="tight")
    print(f"\nWrote {fig_path}")
    print(f"Wrote {outdir/'method_comparison_metrics.csv'}")
    print(f"Wrote {outdir/'method_comparison_master.csv'}")


if __name__ == "__main__":
    main()
