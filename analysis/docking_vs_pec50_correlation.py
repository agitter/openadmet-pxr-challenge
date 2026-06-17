#!/usr/bin/env python3
"""
analysis/04_docking_vs_pec50_correlation.py

Assess how much predictive signal the GNINA ensemble docking scores
(CNNaffinity, minimizedAffinity) carry relative to measured pEC50, using
two independent sets of compounds with known labels:

  Analysis 2: Phase 1 unblinded test compounds (253 compounds)
    - pEC50 from pxr-challenge_TEST_PHASE_1_UNBLINDED.csv
    - docking scores from docking_analysis_extended/docking_receptor_best.csv
      (T### IDs, each compound's best pose across all 62 receptors)

  Analysis 3: Training anchors (89 unique compounds)
    - pEC50 from pxr-challenge_TRAIN.csv
    - docking scores from docking_analysis_extended/docking_receptor_best.csv
      (A#### IDs, each compound's best pose across all 62 receptors)

For each set, computes:
  - Spearman rank correlation (primary, since CNNaffinity/minimizedAffinity
    are not expected to be linearly proportional to pEC50)
  - Pearson correlation (secondary)
  - Kendall tau (matches RAE metric's ranking emphasis)
  - pEC50 range coverage (confirms anchors span the dynamic range)
  - Plots: scatter of best_CNNaffinity vs pEC50, best_minimizedAffinity
    vs pEC50 (saved as PNG)

These correlations set the baseline for what RBFE must improve on, and
confirm the training anchors are well-distributed in score space for
isotonic calibration.

Usage:
    pip install pandas scipy matplotlib
    python analysis/04_docking_vs_pec50_correlation.py \
        --test pxr-challenge_TEST_BLINDED.csv \
        --phase1 pxr-challenge_TEST_PHASE_1_UNBLINDED.csv \
        --train pxr-challenge_TRAIN.csv \
        --receptor-best docking/docking_analysis_extended/docking_receptor_best.csv \
        --outdir docking/docking_analysis_extended
"""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats


def best_scores_per_ligand(receptor_best):
    """From the receptor_best table (one row per cluster x receptor,
    already filtered to best CNNscore pose), get the single best row
    per ligand_id (best CNNscore across all receptors)."""
    idx = receptor_best.groupby("cluster_id")["CNNscore"].idxmax()
    return receptor_best.loc[idx].set_index("cluster_id")


def correlations(x, y, label):
    """Compute and print Spearman, Pearson, Kendall between x and y."""
    mask = ~(np.isnan(x) | np.isnan(y))
    x, y = x[mask], y[mask]
    n = len(x)
    sp = stats.spearmanr(x, y)
    pe = stats.pearsonr(x, y)
    kt = stats.kendalltau(x, y)
    print(f"\n  {label} (n={n}):")
    print(f"    Spearman rho = {sp.statistic:.3f}  (p={sp.pvalue:.2e})")
    print(f"    Pearson  r   = {pe.statistic:.3f}  (p={pe.pvalue:.2e})")
    print(f"    Kendall  tau = {kt.statistic:.3f}  (p={kt.pvalue:.2e})")
    return sp.statistic, pe.statistic, kt.statistic


def scatter_plot(ax, x, y, xlabel, ylabel, title, color):
    ax.scatter(x, y, alpha=0.5, s=20, color=color, edgecolors="none")
    # fit line
    mask = ~(np.isnan(x) | np.isnan(y))
    if mask.sum() > 2:
        m, b = np.polyfit(x[mask], y[mask], 1)
        xr = np.array([x[mask].min(), x[mask].max()])
        ax.plot(xr, m * xr + b, "k--", linewidth=1, alpha=0.7)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    sp = stats.spearmanr(x[mask], y[mask])
    ax.text(0.05, 0.95, f"ρ={sp.statistic:.3f}",
            transform=ax.transAxes, va="top",
            fontsize=9, bbox=dict(boxstyle="round", fc="white", alpha=0.7))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--test", required=True,
                    help="pxr-challenge_TEST_BLINDED.csv")
    ap.add_argument("--phase1", required=True,
                    help="pxr-challenge_TEST_PHASE_1_UNBLINDED.csv")
    ap.add_argument("--train", required=True,
                    help="pxr-challenge_TRAIN.csv")
    ap.add_argument("--receptor-best", required=True,
                    help="docking_analysis_extended/docking_receptor_best.csv")
    ap.add_argument("--outdir", default="docking/docking_analysis_extended")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Load docking scores
    rb = pd.read_csv(args.receptor_best)
    best = best_scores_per_ligand(rb)
    print(f"Loaded docking scores for {len(best)} ligands")

    # ---------------------------------------------------------------
    # Build T### -> Molecule Name mapping (row index of TEST_BLINDED)
    # ---------------------------------------------------------------
    test_df = pd.read_csv(args.test).reset_index(drop=True)
    test_df["ligand_id"] = test_df.index.map(lambda i: f"T{i:03d}")
    test_id_map = test_df.set_index("Molecule Name")["ligand_id"].to_dict()

    # ---------------------------------------------------------------
    # Analysis 2: Phase 1 unblinded test compounds
    # ---------------------------------------------------------------
    print("\n" + "="*60)
    print("ANALYSIS 2: Phase 1 unblinded test compounds (n=253)")
    print("="*60)

    p1 = pd.read_csv(args.phase1)
    p1["ligand_id"] = p1["Molecule Name"].map(test_id_map)
    n_unmapped = p1["ligand_id"].isna().sum()
    if n_unmapped:
        print(f"WARNING: {n_unmapped} Phase 1 compounds not found in test "
              f"set by Molecule Name - check alignment")

    p1_scores = p1.join(best[["CNNscore", "CNNaffinity",
                               "minimizedAffinity", "pdb_id"]],
                         on="ligand_id", how="left")
    n_matched = p1_scores["CNNscore"].notna().sum()
    print(f"Phase 1 compounds with docking scores: {n_matched} / {len(p1)}")

    pec50 = p1_scores["pEC50"].values
    cnn_aff = p1_scores["CNNaffinity"].values
    min_aff = p1_scores["minimizedAffinity"].values
    cnn_score = p1_scores["CNNscore"].values

    print(f"\n  pEC50 range: {np.nanmin(pec50):.2f} - {np.nanmax(pec50):.2f}")
    correlations(cnn_aff, pec50, "CNNaffinity vs pEC50")
    correlations(min_aff, pec50, "minimizedAffinity vs pEC50")
    correlations(cnn_score, pec50, "CNNscore vs pEC50")

    # ---------------------------------------------------------------
    # Analysis 3: Training anchors
    # ---------------------------------------------------------------
    print("\n" + "="*60)
    print("ANALYSIS 3: Training anchors (n=89)")
    print("="*60)

    train_df = pd.read_csv(args.train).reset_index(drop=True)
    train_df["row_idx"] = train_df.index
    train_df["ligand_id"] = train_df["row_idx"].map(lambda i: f"A{i:04d}")

    # only the 89 anchors that were docked
    docked_anchors = best[best.index.str.startswith("A")].copy()
    anchor_ids = docked_anchors.index.tolist()
    print(f"Docked anchor ligand_ids: {len(anchor_ids)}")

    anchor_df = train_df[train_df["ligand_id"].isin(anchor_ids)].copy()
    anchor_df = anchor_df.join(
        docked_anchors[["CNNscore", "CNNaffinity",
                         "minimizedAffinity", "pdb_id"]],
        on="ligand_id", how="left")

    pec50_a = anchor_df["pEC50"].values
    cnn_aff_a = anchor_df["CNNaffinity"].values
    min_aff_a = anchor_df["minimizedAffinity"].values
    cnn_score_a = anchor_df["CNNscore"].values

    print(f"\n  pEC50 range: {np.nanmin(pec50_a):.2f} - {np.nanmax(pec50_a):.2f}")
    print(f"  pEC50 std:   {np.nanstd(pec50_a):.2f}")
    print(f"  CNNaffinity range: {np.nanmin(cnn_aff_a):.2f} - "
          f"{np.nanmax(cnn_aff_a):.2f}")
    correlations(cnn_aff_a, pec50_a, "CNNaffinity vs pEC50")
    correlations(min_aff_a, pec50_a, "minimizedAffinity vs pEC50")
    correlations(cnn_score_a, pec50_a, "CNNscore vs pEC50")

    # ---------------------------------------------------------------
    # Plots
    # ---------------------------------------------------------------
    fig, axes = plt.subplots(2, 2, figsize=(10, 9))
    fig.suptitle("Docking scores vs pEC50", fontsize=13)

    # Phase 1
    scatter_plot(axes[0, 0],
                 cnn_aff, pec50,
                 "CNNaffinity (kcal/mol)", "pEC50",
                 "Phase 1 (n=253): CNNaffinity vs pEC50",
                 color="#2196F3")
    scatter_plot(axes[0, 1],
                 -min_aff, pec50,
                 "-minimizedAffinity (kcal/mol)", "pEC50",
                 "Phase 1 (n=253): -minimizedAffinity vs pEC50",
                 color="#2196F3")

    # Training anchors
    scatter_plot(axes[1, 0],
                 cnn_aff_a, pec50_a,
                 "CNNaffinity (kcal/mol)", "pEC50",
                 "Train anchors (n=89): CNNaffinity vs pEC50",
                 color="#4CAF50")
    scatter_plot(axes[1, 1],
                 -min_aff_a, pec50_a,
                 "-minimizedAffinity (kcal/mol)", "pEC50",
                 "Train anchors (n=89): -minimizedAffinity vs pEC50",
                 color="#4CAF50")

    plt.tight_layout()
    plot_path = outdir / "docking_vs_pec50_correlation.png"
    plt.savefig(plot_path, dpi=150, bbox_inches="tight")
    print(f"\nWrote plot -> {plot_path}")

    # ---------------------------------------------------------------
    # Save merged tables for downstream use
    # ---------------------------------------------------------------
    p1_out = p1_scores[["Molecule Name", "ligand_id", "pEC50",
                          "pEC50_std.error (-log10(molarity))",
                          "CNNscore", "CNNaffinity",
                          "minimizedAffinity", "pdb_id"]].copy()
    p1_out.to_csv(outdir / "phase1_with_docking_scores.csv", index=False)

    anchor_out = anchor_df[["ligand_id", "OCNT_ID", "pEC50",
                              "CNNscore", "CNNaffinity",
                              "minimizedAffinity", "pdb_id"]].copy()
    anchor_out.to_csv(outdir / "train_anchors_with_docking_scores.csv",
                      index=False)

    print(f"Wrote phase1_with_docking_scores.csv "
          f"({len(p1_out)} rows)")
    print(f"Wrote train_anchors_with_docking_scores.csv "
          f"({len(anchor_out)} rows)")


if __name__ == "__main__":
    main()
