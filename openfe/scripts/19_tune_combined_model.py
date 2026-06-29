#!/usr/bin/env python3
"""
openfe/scripts/19_tune_combined_model.py

Tune and select a combined docking+RBFE model for pEC50 prediction,
validated on the unblinded Phase 1 set using the official competition
RAE metric (bootstrapped, seed 0). After selection, the WINNING model is
refit on all 253 Phase 1 compounds and its COMPLETE parameters are
serialized to chosen_model.json, so the Set2 application script needs no
fitting logic - it just loads and applies the frozen transform.

Design:
  - Each model is defined by a fit_* function (returns a params dict) and
    an apply_* function (params + features -> predictions). The SAME
    functions are used inside CV folds and for the final full-data refit,
    so the deployed model matches exactly what was cross-validated.
  - Evaluate all models on all 253 Phase 1 compounds via shared 5-fold CV.
  - w=0 when RBFE absent (pure docking); rbfe_only uses train-mean
    placeholder where RBFE missing.
  - Pool out-of-fold predictions, bootstrap RAE 1000x (seed 0).
  - Winner refit on all 253; full params saved.

Gate features: path_ddg_error, min_overlap_on_path, max_edge_error_on_path,
               n_hops, std_CNNscore_across_receptors

Usage:
    python openfe/scripts/19_tune_combined_model.py \
        --rbfe-pred openfe/rbfe_predictions.csv \
        --receptor-best docking/docking_analysis_extended/docking_receptor_best.csv \
        --cluster-summary docking/docking_analysis_extended/docking_cluster_summary.csv \
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
from scipy.optimize import minimize
from scipy.stats import spearmanr


def rae(y_true, y_pred):
    """Official competition RAE (from challenge config.py)."""
    return np.sum(np.abs(y_true - y_pred)) / np.sum(
        np.abs(y_true - np.mean(y_true)))


BOOTSTRAP_SEED = 0
N_BOOTSTRAP = 1000
GATE_FEATURES = ["path_ddg_error", "min_overlap_on_path",
                 "max_edge_error_on_path", "n_hops",
                 "std_CNNscore_across_receptors"]
MODELS = ["docking_only", "rbfe_only", "fixed_blend",
          "threshold_gate", "sigmoid_gate"]


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
    return float(vals.mean()), float(vals.std()), float(rae(y_true, y_pred))


def linfit(x, y):
    A = np.vstack([x, np.ones_like(x)]).T
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    return float(coef[0]), float(coef[1])


def sigmoid(z):
    return 1.0 / (1.0 + np.exp(-np.clip(z, -30, 30)))


# ---------------------------------------------------------------------------
# Fit / apply functions. Each fit_* returns a JSON-serializable params dict;
# each apply_* takes (params, frame) and returns predictions. The frame must
# carry: CNNaffinity, pred_pEC50_raw, has_rbfe, and the GATE_FEATURES.
# Calibrations are shared and passed in via `cal`.
# ---------------------------------------------------------------------------

def fit_calibrations(train):
    a_d, b_d = linfit(train["CNNaffinity"].values, train["exp_pEC50"].values)
    tr_rb = train[train["has_rbfe"]]
    if len(tr_rb) >= 5:
        a_r, b_r = linfit(tr_rb["pred_pEC50_raw"].values,
                          tr_rb["exp_pEC50"].values)
    else:
        a_r, b_r = 0.0, float(train["exp_pEC50"].mean())
    return {"dock_slope": a_d, "dock_intercept": b_d,
            "rbfe_slope": a_r, "rbfe_intercept": b_r,
            "train_mean": float(train["exp_pEC50"].mean())}


def dock_cal(cal, frame):
    return cal["dock_slope"] * frame["CNNaffinity"].values + cal["dock_intercept"]


def rbfe_cal(cal, frame):
    raw = frame["pred_pEC50_raw"].values
    out = cal["rbfe_slope"] * raw + cal["rbfe_intercept"]
    return np.where(frame["has_rbfe"].values, out, cal["train_mean"])


def feat_matrix(frame):
    return np.nan_to_num(frame[GATE_FEATURES].values.astype(float), nan=0.0)


# ---- docking_only ----
def fit_docking_only(train, cal):
    return {}

def apply_docking_only(params, cal, frame):
    return dock_cal(cal, frame)


# ---- rbfe_only ----
def fit_rbfe_only(train, cal):
    return {}

def apply_rbfe_only(params, cal, frame):
    return rbfe_cal(cal, frame)


# ---- fixed_blend ----
def fit_fixed_blend(train, cal):
    d = dock_cal(cal, train)
    r = rbfe_cal(cal, train)
    has = train["has_rbfe"].values
    t = train["exp_pEC50"].values
    best_w, best = 0.0, np.inf
    for w in np.linspace(0, 1, 21):
        wv = np.where(has, w, 0.0)
        cur = rae(t, wv * r + (1 - wv) * d)
        if cur < best:
            best, best_w = cur, w
    return {"w": float(best_w)}

def apply_fixed_blend(params, cal, frame):
    d = dock_cal(cal, frame); r = rbfe_cal(cal, frame)
    wv = np.where(frame["has_rbfe"].values, params["w"], 0.0)
    return wv * r + (1 - wv) * d


# ---- threshold_gate ----
def fit_threshold_gate(train, cal):
    d = dock_cal(cal, train); r = rbfe_cal(cal, train)
    has = train["has_rbfe"].values
    t = train["exp_pEC50"].values
    err = train["path_ddg_error"].values
    err = np.nan_to_num(err, nan=99.0)
    best, best_cfg = np.inf, (1.5, 3.0, 0.0, 0.5)
    for t_low in [1.0, 1.5, 2.0]:
        for t_high in [2.5, 3.0, 3.5]:
            for w_lo in [0.0, 0.1, 0.2]:
                for w_hi in [0.3, 0.5, 0.7]:
                    w = np.where(err < t_low, w_hi,
                                 np.where(err < t_high, (w_hi+w_lo)/2, w_lo))
                    w = np.where(has, w, 0.0)
                    cur = rae(t, w*r + (1-w)*d)
                    if cur < best:
                        best, best_cfg = cur, (t_low, t_high, w_lo, w_hi)
    return {"t_low": best_cfg[0], "t_high": best_cfg[1],
            "w_lo": best_cfg[2], "w_hi": best_cfg[3]}

def apply_threshold_gate(params, cal, frame):
    d = dock_cal(cal, frame); r = rbfe_cal(cal, frame)
    err = np.nan_to_num(frame["path_ddg_error"].values, nan=99.0)
    tl, th, wl, wh = (params["t_low"], params["t_high"],
                      params["w_lo"], params["w_hi"])
    w = np.where(err < tl, wh, np.where(err < th, (wh+wl)/2, wl))
    w = np.where(frame["has_rbfe"].values, w, 0.0)
    return w * r + (1 - w) * d


# ---- sigmoid_gate ----
def fit_sigmoid_gate(train, cal, rng):
    d = dock_cal(cal, train); r = rbfe_cal(cal, train)
    has = train["has_rbfe"].values
    t = train["exp_pEC50"].values
    raw = train[GATE_FEATURES].values.astype(float)
    mu = np.nanmean(raw, axis=0); sd = np.nanstd(raw, axis=0)
    sd[sd == 0] = 1.0
    F = np.nan_to_num((np.nan_to_num(raw, nan=0.0) - mu) / sd, nan=0.0)

    def obj(beta):
        w = sigmoid(F @ beta[1:] + beta[0])
        w = np.where(has, w, 0.0)
        return rae(t, w*r + (1-w)*d) + 0.01 * np.sum(beta[1:]**2)

    best_beta, best = None, np.inf
    for _ in range(5):
        b0 = rng.normal(0, 0.5, size=len(GATE_FEATURES)+1)
        res = minimize(obj, b0, method="Nelder-Mead",
                       options={"maxiter": 2000, "xatol": 1e-4, "fatol": 1e-4})
        if res.fun < best:
            best, best_beta = res.fun, res.x
    return {"beta": best_beta.tolist(), "mu": mu.tolist(), "sd": sd.tolist()}

def apply_sigmoid_gate(params, cal, frame):
    d = dock_cal(cal, frame); r = rbfe_cal(cal, frame)
    raw = frame[GATE_FEATURES].values.astype(float)
    mu = np.array(params["mu"]); sd = np.array(params["sd"])
    F = np.nan_to_num((np.nan_to_num(raw, nan=0.0) - mu) / sd, nan=0.0)
    beta = np.array(params["beta"])
    w = sigmoid(F @ beta[1:] + beta[0])
    w = np.where(frame["has_rbfe"].values, w, 0.0)
    return w * r + (1 - w) * d


FIT = {"docking_only": fit_docking_only, "rbfe_only": fit_rbfe_only,
       "fixed_blend": fit_fixed_blend, "threshold_gate": fit_threshold_gate,
       "sigmoid_gate": fit_sigmoid_gate}
APPLY = {"docking_only": apply_docking_only, "rbfe_only": apply_rbfe_only,
         "fixed_blend": apply_fixed_blend, "threshold_gate": apply_threshold_gate,
         "sigmoid_gate": apply_sigmoid_gate}


def fit_model(name, train, cal, rng):
    if name == "sigmoid_gate":
        return FIT[name](train, cal, rng)
    return FIT[name](train, cal)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rbfe-pred", default="openfe/rbfe_predictions.csv")
    ap.add_argument("--receptor-best",
                    default="docking/docking_analysis_extended/"
                            "docking_receptor_best.csv")
    ap.add_argument("--cluster-summary",
                    default="docking/docking_analysis_extended/"
                            "docking_cluster_summary.csv")
    ap.add_argument("--phase1",
                    default="data/pxr-challenge_TEST_PHASE_1_UNBLINDED.csv")
    ap.add_argument("--n-folds", type=int, default=5)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--outdir", default="openfe")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    rbfe = pd.read_csv(args.rbfe_pred)
    best = pd.read_csv(args.receptor_best)
    csum = pd.read_csv(args.cluster_summary)
    phase1 = pd.read_csv(args.phase1)

    dock = (best.sort_values("CNNaffinity", ascending=False)
                .groupby("ligand_name", as_index=False).first()
                [["ligand_name", "CNNaffinity"]]
                .rename(columns={"ligand_name": "Molecule Name"}))
    if "std_CNNscore_across_receptors" in csum.columns:
        cstd = csum[["ligand_name", "std_CNNscore_across_receptors"]].rename(
            columns={"ligand_name": "Molecule Name"})
    else:
        cstd = pd.DataFrame(columns=["Molecule Name",
                                     "std_CNNscore_across_receptors"])

    df = phase1[["Molecule Name", "pEC50"]].rename(
        columns={"pEC50": "exp_pEC50"}).copy()
    df = df.merge(dock, on="Molecule Name", how="left")
    df = df.merge(cstd, on="Molecule Name", how="left")
    rcols = ["Molecule Name", "pred_pEC50_raw", "path_ddg_error",
             "min_overlap_on_path", "max_edge_error_on_path", "n_hops"]
    df = df.merge(rbfe[rcols], on="Molecule Name", how="left")
    df["has_rbfe"] = df["pred_pEC50_raw"].notna()
    df = df[df["CNNaffinity"].notna()].reset_index(drop=True)
    print(f"Phase 1 compounds: {len(df)}  "
          f"(with RBFE: {int(df['has_rbfe'].sum())}, "
          f"docking-only: {int((~df['has_rbfe']).sum())})")

    rng = np.random.default_rng(args.seed)
    order = rng.permutation(len(df))
    folds = np.array_split(order, args.n_folds)

    oof = {m: np.full(len(df), np.nan) for m in MODELS}
    for fold_id in range(args.n_folds):
        test_idx = folds[fold_id]
        train_idx = np.concatenate(
            [folds[j] for j in range(args.n_folds) if j != fold_id])
        tr, te = df.iloc[train_idx], df.iloc[test_idx]
        cal = fit_calibrations(tr)
        for m in MODELS:
            params = fit_model(m, tr, cal, rng)
            oof[m][test_idx] = APPLY[m](params, cal, te)

    y_true = df["exp_pEC50"].values
    rows = []
    for m in MODELS:
        mean_rae, std_rae, point_rae = bootstrap_rae(y_true, oof[m])
        rho, _ = spearmanr(y_true, oof[m])
        mae = float(np.mean(np.abs(y_true - oof[m])))
        rows.append({"model": m, "rae_boot_mean": mean_rae,
                     "rae_boot_std": std_rae, "rae_point": point_rae,
                     "spearman": rho, "mae": mae})
    metrics = pd.DataFrame(rows).sort_values("rae_boot_mean")
    metrics.to_csv(outdir / "model_comparison_cv.csv", index=False)
    print("\n" + "=" * 64)
    print("MODEL COMPARISON (5-fold CV, bootstrapped RAE, lower=better)")
    print("=" * 64)
    print(metrics.to_string(index=False))

    winner = metrics.iloc[0]["model"]
    print(f"\nBest model: {winner} "
          f"(RAE {metrics.iloc[0]['rae_boot_mean']:.3f} "
          f"+/- {metrics.iloc[0]['rae_boot_std']:.3f})")

    # ---- Visualizations ----
    fig, axes = plt.subplots(2, 3, figsize=(17, 10))
    ax = axes[0, 0]
    ax.barh(metrics["model"], metrics["rae_boot_mean"],
            xerr=metrics["rae_boot_std"], color="#6baed6",
            edgecolor="k", capsize=3)
    ax.axvline(1.0, color="r", ls="--", lw=0.8, label="mean-predictor")
    ax.set_xlabel("Bootstrapped RAE (lower better)")
    ax.set_title("Model comparison (Phase 1 CV)")
    ax.legend(fontsize=8); ax.invert_yaxis()

    positions = [(0, 1), (0, 2), (1, 0), (1, 1), (1, 2)]
    covered = df["has_rbfe"].values
    for (m, pos) in zip(MODELS, positions):
        ax = axes[pos]; yp = oof[m]
        ax.scatter(yp[~covered], y_true[~covered], s=14, alpha=0.5,
                   c="#74c476", label="docking-only")
        ax.scatter(yp[covered], y_true[covered], s=14, alpha=0.5,
                   c="#fd8d3c", label="has RBFE")
        lims = [min(y_true.min(), np.nanmin(yp))-0.5,
                max(y_true.max(), np.nanmax(yp))+0.5]
        ax.plot(lims, lims, "k--", lw=0.7, alpha=0.6)
        ax.set_xlim(lims); ax.set_ylim(lims)
        ax.set_xlabel(f"{m} predicted pEC50")
        ax.set_ylabel("experimental pEC50")
        mrow = metrics[metrics["model"] == m].iloc[0]
        ax.set_title(f"{m}\nRAE={mrow['rae_boot_mean']:.3f} "
                     f"rho={mrow['spearman']:.3f}")
        ax.legend(fontsize=7)
    plt.tight_layout()
    fig_path = outdir / "model_tuning_diagnostics.png"
    plt.savefig(fig_path, dpi=130, bbox_inches="tight")
    print(f"\nWrote {fig_path}")

    # ---- Refit winner on ALL 253 and serialize COMPLETE model ----
    cal_full = fit_calibrations(df)
    winner_params = fit_model(winner, df, cal_full, rng)
    model_json = {
        "winner": winner,
        "gate_features": GATE_FEATURES,
        "calibrations": cal_full,
        "winner_params": winner_params,
        "docking_input": "best_CNNaffinity_per_compound",
        "rbfe_input": "pred_pEC50_raw",
        "rae_cv_mean": float(metrics.iloc[0]["rae_boot_mean"]),
        "rae_cv_std": float(metrics.iloc[0]["rae_boot_std"]),
        "apply_rule": (
            "pred = dock_cal where no RBFE; otherwise apply winner model. "
            "dock_cal = dock_slope*CNNaffinity + dock_intercept. "
            "rbfe_cal = rbfe_slope*pred_pEC50_raw + rbfe_intercept."),
    }
    with open(outdir / "chosen_model.json", "w") as f:
        json.dump(model_json, f, indent=2)
    print(f"Wrote {outdir/'chosen_model.json'} "
          f"(complete frozen model: {winner})")


if __name__ == "__main__":
    main()
