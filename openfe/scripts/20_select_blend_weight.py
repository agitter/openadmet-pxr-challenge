#!/usr/bin/env python3
"""
openfe/scripts/20_select_blend_weight.py

Targeted selection of the docking/RBFE blend weight w for the final
submission, chosen to give RBFE a meaningful voice (for methods
evaluation) while keeping RAE near optimal.

The blend (where RBFE is present): pred = w*rbfe_cal + (1-w)*dock_cal.
Where RBFE is absent: pred = dock_cal (w effectively 0).

Sweeps w on all 253 Phase 1 compounds, computes bootstrapped RAE
(mean +/- std, seed 0) at each w, and marks three candidate operating
points:
  1. RAE-optimal w (lowest mean RAE)
  2. Largest w within RAE tolerance (default +0.005) of the optimum
  3. "Elbow": where the RAE-vs-w slope starts climbing steeply

Saves the chosen model (calibrations + selected w) to chosen_model.json
in the same frozen format script 19 used, ready for the Set2 applier.

Usage:
    python openfe/scripts/20_select_blend_weight.py \
        --rbfe-pred openfe/rbfe_predictions.csv \
        --receptor-best docking/docking_analysis_extended/docking_receptor_best.csv \
        --phase1 data/pxr-challenge_TEST_PHASE_1_UNBLINDED.csv \
        --choice elbow \
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

BOOTSTRAP_SEED = 0
N_BOOTSTRAP = 1000


def rae(y_true, y_pred):
    return np.sum(np.abs(y_true - y_pred)) / np.sum(
        np.abs(y_true - np.mean(y_true)))


def bootstrap_rae(y_true, y_pred, n=N_BOOTSTRAP, seed=BOOTSTRAP_SEED):
    rng = np.random.default_rng(seed=seed)
    m = len(y_true)
    idx = rng.choice(m, size=(n, m), replace=True)
    vals = []
    for ix in idx:
        yt, yp = y_true[ix], y_pred[ix]
        denom = np.sum(np.abs(yt - np.mean(yt)))
        if denom > 0:
            vals.append(np.sum(np.abs(yt - yp)) / denom)
    vals = np.array(vals)
    return vals.mean(), vals.std()


def linfit(x, y):
    A = np.vstack([x, np.ones_like(x)]).T
    c, *_ = np.linalg.lstsq(A, y, rcond=None)
    return float(c[0]), float(c[1])


def find_elbow(ws, raes):
    """Elbow via max distance from the chord connecting the curve's
    endpoints (Kneedle-style geometric method)."""
    x = np.asarray(ws); y = np.asarray(raes)
    x0, x1 = x[0], x[-1]; y0, y1 = y[0], y[-1]
    # distance of each point from the line through (x0,y0)-(x1,y1)
    num = np.abs((y1 - y0) * x - (x1 - x0) * y + x1 * y0 - y1 * x0)
    den = np.hypot(y1 - y0, x1 - x0)
    dist = num / den
    # only consider the rising portion after the minimum
    imin = int(np.argmin(y))
    if imin >= len(x) - 2:
        return x[imin]
    dist[:imin] = -np.inf
    return float(x[int(np.argmax(dist))])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rbfe-pred", default="openfe/rbfe_predictions.csv")
    ap.add_argument("--receptor-best",
                    default="docking/docking_analysis_extended/"
                            "docking_receptor_best.csv")
    ap.add_argument("--phase1",
                    default="data/pxr-challenge_TEST_PHASE_1_UNBLINDED.csv")
    ap.add_argument("--tolerance", type=float, default=0.005,
                    help="RAE tolerance above optimum for candidate 2")
    ap.add_argument("--choice", default="elbow",
                    choices=["optimal", "tolerance", "elbow"],
                    help="Which candidate w to freeze into chosen_model.json")
    ap.add_argument("--outdir", default="openfe")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    rbfe = pd.read_csv(args.rbfe_pred)
    best = pd.read_csv(args.receptor_best)
    phase1 = pd.read_csv(args.phase1)

    dock = (best.sort_values("CNNaffinity", ascending=False)
                .groupby("ligand_name", as_index=False).first()
                [["ligand_name", "CNNaffinity"]]
                .rename(columns={"ligand_name": "Molecule Name"}))
    df = phase1[["Molecule Name", "pEC50"]].rename(
        columns={"pEC50": "exp_pEC50"})
    df = df.merge(dock, on="Molecule Name", how="left")
    df = df.merge(rbfe[["Molecule Name", "pred_pEC50_raw"]],
                  on="Molecule Name", how="left")
    df["has_rbfe"] = df["pred_pEC50_raw"].notna()
    df = df[df["CNNaffinity"].notna()].reset_index(drop=True)
    print(f"Phase 1: {len(df)} compounds "
          f"({int(df['has_rbfe'].sum())} with RBFE)")

    # Calibrations on all 253
    a_d, b_d = linfit(df["CNNaffinity"].values, df["exp_pEC50"].values)
    rb = df[df["has_rbfe"]]
    a_r, b_r = linfit(rb["pred_pEC50_raw"].values, rb["exp_pEC50"].values)
    train_mean = float(df["exp_pEC50"].mean())

    d = a_d * df["CNNaffinity"].values + b_d
    r = np.where(df["has_rbfe"].values,
                 a_r * df["pred_pEC50_raw"].values + b_r, train_mean)
    has = df["has_rbfe"].values
    t = df["exp_pEC50"].values

    ws = np.linspace(0, 1, 41)
    means, stds = [], []
    for w in ws:
        wv = np.where(has, w, 0.0)
        pred = wv * r + (1 - wv) * d
        mean_rae, std_rae = bootstrap_rae(t, pred)
        means.append(mean_rae); stds.append(std_rae)
    means = np.array(means); stds = np.array(stds)

    # Candidates
    w_opt = float(ws[int(np.argmin(means))])
    opt_rae = means.min()
    within = ws[means <= opt_rae + args.tolerance]
    w_tol = float(within.max()) if len(within) else w_opt
    w_elbow = find_elbow(ws, means)

    candidates = {"optimal": w_opt, "tolerance": w_tol, "elbow": w_elbow}
    print("\nCandidate blend weights:")
    for name, w in candidates.items():
        wv = np.where(has, w, 0.0)
        mr, sr = bootstrap_rae(t, wv * r + (1 - wv) * d)
        print(f"  {name:10s}: w={w:.3f}  RAE={mr:.4f} +/- {sr:.4f}")

    chosen_w = candidates[args.choice]
    print(f"\nChosen ({args.choice}): w={chosen_w:.3f}")

    # ---- Plot ----
    fig, ax = plt.subplots(figsize=(9, 6))
    ax.plot(ws, means, "-", color="#2171b5", lw=1.5)
    ax.fill_between(ws, means - stds, means + stds, alpha=0.2,
                    color="#6baed6", label="bootstrap +/- 1 SD")
    ax.axhline(1.0, color="r", ls="--", lw=0.8, label="mean-predictor")
    colors = {"optimal": "#31a354", "tolerance": "#e6550d",
              "elbow": "#756bb1"}
    for name, w in candidates.items():
        ax.axvline(w, color=colors[name], ls=":", lw=1.5,
                   label=f"{name} (w={w:.2f})")
    ax.scatter([chosen_w], [means[np.argmin(np.abs(ws - chosen_w))]],
               s=120, edgecolor="k", facecolor="gold", zorder=5,
               label=f"CHOSEN: {args.choice}")
    ax.set_xlabel("Blend weight w (on RBFE, where present)")
    ax.set_ylabel("Bootstrapped RAE")
    ax.set_title("Docking/RBFE blend weight selection (Phase 1)\n"
                 "w=0 pure docking; w=1 pure RBFE where available")
    ax.legend(fontsize=8)
    fig_path = outdir / "blend_weight_selection.png"
    plt.savefig(fig_path, dpi=130, bbox_inches="tight")
    print(f"Wrote {fig_path}")

    # ---- Save chosen model (frozen, same format as script 19) ----
    model_json = {
        "winner": "fixed_blend",
        "selection_method": args.choice,
        "calibrations": {
            "dock_slope": a_d, "dock_intercept": b_d,
            "rbfe_slope": a_r, "rbfe_intercept": b_r,
            "train_mean": train_mean},
        "winner_params": {"w": chosen_w},
        "candidate_weights": candidates,
        "docking_input": "best_CNNaffinity_per_compound",
        "rbfe_input": "pred_pEC50_raw",
        "apply_rule": (
            "pred = w*rbfe_cal + (1-w)*dock_cal where RBFE present, "
            "else dock_cal. dock_cal=dock_slope*CNNaffinity+dock_intercept; "
            "rbfe_cal=rbfe_slope*pred_pEC50_raw+rbfe_intercept."),
        "note": (
            "Blend weight deliberately set above the RAE optimum to give "
            "RBFE a meaningful voice for methods evaluation. RBFE is "
            "expected to underperform docking on this target; see writeup."),
    }
    with open(outdir / "chosen_model.json", "w") as f:
        json.dump(model_json, f, indent=2)
    print(f"Wrote {outdir/'chosen_model.json'} (fixed_blend, w={chosen_w:.3f})")


if __name__ == "__main__":
    main()
