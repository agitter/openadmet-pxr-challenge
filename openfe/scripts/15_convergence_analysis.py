#!/usr/bin/env python3
"""
openfe/scripts/15_convergence_analysis.py

Campaign-wide convergence analysis. Extracts MBAR statistical error,
lambda-state overlap, forward/reverse agreement, and replica-exchange
mixing for EVERY completed leg, then summarizes convergence quality
and correlates it with |ddG| to identify unreliable edges.

The convergence diagnostics live in the Analysis unit's outputs inside
result.json:
  - unit_estimate, unit_estimate_error  (dG and MBAR error)
  - unit_mbar_overlap                   (overlap matrix; min adjacent
                                         overlap is the key scalar)
  - forward_and_reverse_energies        (convergence over time)
  - replica_exchange_statistics         (mixing)
  - production_iterations               (sampling completed)

Outputs:
  openfe/leg_convergence.csv  - per-leg diagnostics
  openfe/edge_convergence.csv - per-edge, joined with ddG

Usage:
    python openfe/scripts/15_convergence_analysis.py \
        --production-dir openfe/production \
        --edge-results openfe/all_edge_results.csv \
        --outdir openfe
"""

import argparse
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd


def get_analysis_outputs(result_path):
    """Return the Analysis unit's outputs dict, or None."""
    if not result_path.exists() or result_path.stat().st_size <= 100:
        return None
    try:
        d = json.loads(result_path.read_text())
    except Exception:
        return None
    # Look in unit_results for the Analysis unit
    for k, v in d.get("unit_results", {}).items():
        name = v.get("name", "")
        if "Analysis" in name:
            return v.get("outputs", {})
    # Also check protocol_result.data list
    pr_data = d.get("protocol_result", {}).get("data", {})
    for hashkey, val in pr_data.items():
        if isinstance(val, list):
            for item in val:
                if isinstance(item, dict) and "Analysis" in item.get("name", ""):
                    return item.get("outputs", {})
    return None


def scalar(x):
    """Extract a float magnitude from an openfe quantity dict or number."""
    if x is None:
        return None
    if isinstance(x, dict):
        return x.get("magnitude")
    if isinstance(x, (int, float)):
        return float(x)
    return None


def overlap_scalar(overlap):
    """The mbar overlap is a dict with a 'scalar' field - the overlap
    measure (higher = better phase-space overlap between lambda states).
    Values near 0 indicate poor overlap / unreliable estimate."""
    if overlap is None:
        return None
    if isinstance(overlap, dict):
        return overlap.get("scalar")
    return None


def extract_leg(result_path):
    """Pull convergence diagnostics for one leg."""
    out = get_analysis_outputs(result_path)
    if out is None:
        return None
    diag = {
        "dG": scalar(out.get("unit_estimate")),
        "dG_error": scalar(out.get("unit_estimate_error")),
        "production_iterations": out.get("production_iterations"),
        "equilibration_iterations": out.get("equilibration_iterations"),
        "min_overlap": overlap_scalar(out.get("unit_mbar_overlap")),
    }
    # Forward/reverse agreement: last forward vs last reverse estimate
    fr = out.get("forward_and_reverse_energies")
    diag["fwd_rev_diff"] = None
    if isinstance(fr, dict):
        try:
            fwd = fr.get("forward")
            rev = fr.get("reverse")
            if fwd and rev:
                fwd_last = scalar(fwd[-1]) if isinstance(fwd, list) else scalar(fwd)
                rev_last = scalar(rev[-1]) if isinstance(rev, list) else scalar(rev)
                if fwd_last is not None and rev_last is not None:
                    diag["fwd_rev_diff"] = abs(fwd_last - rev_last)
        except Exception:
            pass
    return diag


def parse_edge_ligands(edge_name):
    body = edge_name[len("rbfe_"):] if edge_name.startswith("rbfe_") else edge_name
    parts = body.split("_")
    if len(parts) == 2:
        return parts[0], parts[1]
    return None, None


def find_leg_result(prod_dir, cluster_id, ligA, ligB, leg):
    cluster_dir = prod_dir / str(cluster_id)
    if not cluster_dir.exists():
        return None
    for a, b in [(ligA, ligB), (ligB, ligA)]:
        name = f"rbfe_{a}_{leg}_{b}_{leg}"
        d = cluster_dir / name
        if d.exists():
            for sub in ["result.json", "quickrun_output/result.json",
                        "salvage/result.json"]:
                rp = d / sub
                if rp.exists() and rp.stat().st_size > 100:
                    # only return if it has a valid estimate
                    try:
                        data = json.loads(rp.read_text())
                        est = data.get("estimate")
                        if est and est.get("magnitude") is not None:
                            return rp
                    except Exception:
                        pass
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--production-dir", default="openfe/production")
    ap.add_argument("--edge-results", default="openfe/all_edge_results.csv")
    ap.add_argument("--outdir", default="openfe")
    args = ap.parse_args()

    prod_dir = Path(args.production_dir)
    outdir = Path(args.outdir)
    edge_df = pd.read_csv(args.edge_results)
    done = edge_df[edge_df["both_done"]].copy()

    edge_rows = []
    leg_rows = []

    for _, row in done.iterrows():
        cluster_id = row["cluster_id"]
        edge = row["edge"]
        ligA, ligB = parse_edge_ligands(edge)
        rec = {"cluster_id": cluster_id, "edge": edge, "ddg": row["ddg"]}

        for leg in ["complex", "solvent"]:
            rp = find_leg_result(prod_dir, cluster_id, ligA, ligB, leg)
            if rp is None:
                continue
            diag = extract_leg(rp)
            if diag is None:
                continue
            leg_rows.append({"cluster_id": cluster_id, "edge": edge,
                             "leg": leg, **diag})
            for key, val in diag.items():
                rec[f"{leg}_{key}"] = val

        edge_rows.append(rec)

    leg_df = pd.DataFrame(leg_rows)
    edge_conv = pd.DataFrame(edge_rows)
    leg_df.to_csv(outdir / "leg_convergence.csv", index=False)
    edge_conv.to_csv(outdir / "edge_convergence.csv", index=False)

    print("=" * 60)
    print("CONVERGENCE ANALYSIS (all completed legs)")
    print("=" * 60)
    print(f"Legs with diagnostics: {len(leg_df)}")

    if len(leg_df):
        print(f"\nMBAR dG_error (kcal/mol):")
        print(leg_df["dG_error"].describe().to_string())
        print(f"\nMin adjacent lambda overlap:")
        print(leg_df["min_overlap"].describe().to_string())

        # Correlate convergence with |ddG|
        edge_conv["absddg"] = edge_conv["ddg"].abs()
        # Worst error across the two legs
        err_cols = [c for c in edge_conv.columns if c.endswith("_dG_error")]
        ovl_cols = [c for c in edge_conv.columns if c.endswith("_min_overlap")]
        if err_cols:
            edge_conv["max_leg_error"] = edge_conv[err_cols].max(axis=1)
        if ovl_cols:
            edge_conv["min_leg_overlap"] = edge_conv[ovl_cols].min(axis=1)

        print("\n" + "=" * 60)
        print("CONVERGENCE vs |ddG|")
        print("=" * 60)
        for thresh in [3, 5, 7]:
            big = edge_conv[edge_conv["absddg"] > thresh]
            small = edge_conv[edge_conv["absddg"] <= thresh]
            if "max_leg_error" in edge_conv.columns:
                print(f"\n|ddG| > {thresh} (n={len(big)}): "
                      f"mean MBAR err={big['max_leg_error'].mean():.3f}")
                print(f"|ddG| <= {thresh} (n={len(small)}): "
                      f"mean MBAR err={small['max_leg_error'].mean():.3f}")
            if "min_leg_overlap" in edge_conv.columns:
                print(f"  |ddG|>{thresh} mean min-overlap="
                      f"{big['min_leg_overlap'].mean():.4f}  "
                      f"|ddG|<={thresh} mean min-overlap="
                      f"{small['min_leg_overlap'].mean():.4f}")

    print(f"\nWrote {outdir/'leg_convergence.csv'} and "
          f"{outdir/'edge_convergence.csv'}")


if __name__ == "__main__":
    main()
