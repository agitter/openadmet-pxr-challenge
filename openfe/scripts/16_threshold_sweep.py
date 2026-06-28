#!/usr/bin/env python3
"""
openfe/scripts/16_threshold_sweep.py

Use the unblinded Phase 1 compounds as ground truth to choose an
edge-reliability threshold empirically.

For each candidate ddG_error threshold:
  1. Keep only edges whose propagated MBAR ddG_error <= threshold
  2. Rebuild per-cluster graphs and re-propagate pEC50 from anchors
  3. Compare RBFE predictions vs known Phase 1 pEC50 (Spearman, Pearson, MAE)
  4. Report how many compounds remain connected (coverage)

This trades coverage against accuracy and lets the data pick the
operating point that maximizes agreement with Phase 1 while preserving
as many predictions as possible.

Inputs:
  - openfe/edge_convergence.csv  (per-edge ddg + complex/solvent dG_error)
  - openfe/test_full_with_clusters_and_anchors.csv
  - data/pxr-challenge_TRAIN.csv
  - data/pxr-challenge_TEST_PHASE_1_UNBLINDED.csv

Output:
  - openfe/threshold_sweep.csv

Usage:
    python openfe/scripts/16_threshold_sweep.py \
        --edge-convergence openfe/edge_convergence.csv \
        --test-full openfe/test_full_with_clusters_and_anchors.csv \
        --train data/pxr-challenge_TRAIN.csv \
        --phase1 data/pxr-challenge_TEST_PHASE_1_UNBLINDED.csv \
        --outdir openfe
"""

import argparse
from collections import defaultdict, deque
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

RT_LN10 = 1.3638


def parse_edge_ligands(edge_name):
    body = edge_name[len("rbfe_"):] if edge_name.startswith("rbfe_") else edge_name
    parts = body.split("_")
    if len(parts) == 2:
        return parts[0], parts[1]
    return None, None


def build_graph(edges):
    """edges: list of (a, b, ddg). Directed-ddG adjacency."""
    graph = defaultdict(dict)
    for a, b, ddg in edges:
        graph[a][b] = ddg
        graph[b][a] = -ddg
    return graph


def bfs_paths(graph, start):
    visited = {start: (0.0, 0)}
    queue = deque([start])
    while queue:
        node = queue.popleft()
        curr_ddg, curr_hops = visited[node]
        for neighbor, ddg in graph[node].items():
            if neighbor not in visited:
                visited[neighbor] = (curr_ddg + ddg, curr_hops + 1)
                queue.append(neighbor)
    return visited


def propagate(edge_df, anchor_pec50, max_error=None):
    """Propagate pEC50 from anchors, optionally filtering edges by
    ddg_error <= max_error. Returns dict compound -> pred_pEC50."""
    preds = {}
    for cluster_id, grp in edge_df.groupby("cluster_id"):
        edges = []
        for _, row in grp.iterrows():
            if max_error is not None and row["ddg_error"] > max_error:
                continue
            a, b = parse_edge_ligands(row["edge"])
            if a is None:
                continue
            edges.append((a, b, row["ddg"]))
        if not edges:
            continue
        graph = build_graph(edges)
        nodes = set(graph.keys())
        anchors = {n for n in nodes if n in anchor_pec50}
        test_compounds = {n for n in nodes if n.startswith("OADMET")}
        if not anchors:
            continue
        best = {}
        for anchor in anchors:
            paths = bfs_paths(graph, anchor)
            for comp in test_compounds:
                if comp in paths:
                    path_ddg, hops = paths[comp]
                    if comp not in best or hops < best[comp][1]:
                        best[comp] = (path_ddg, hops, anchor)
        for comp, (path_ddg, hops, anchor) in best.items():
            preds[comp] = anchor_pec50[anchor] - path_ddg / RT_LN10
    return preds


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--edge-convergence", default="openfe/edge_convergence.csv")
    ap.add_argument("--test-full",
                    default="openfe/test_full_with_clusters_and_anchors.csv")
    ap.add_argument("--train", default="data/pxr-challenge_TRAIN.csv")
    ap.add_argument("--phase1",
                    default="data/pxr-challenge_TEST_PHASE_1_UNBLINDED.csv")
    ap.add_argument("--outdir", default="openfe")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    edge_df = pd.read_csv(args.edge_convergence)
    test_df = pd.read_csv(args.test_full)
    train_df = pd.read_csv(args.train)
    phase1_df = pd.read_csv(args.phase1)

    # Compute propagated ddg_error per edge from the two leg errors
    ecols = [c for c in edge_df.columns if c.endswith("_dG_error")]
    if len(ecols) >= 2:
        edge_df["ddg_error"] = np.sqrt(
            edge_df[ecols[0]].fillna(0)**2 + edge_df[ecols[1]].fillna(0)**2)
    else:
        edge_df["ddg_error"] = edge_df[ecols[0]].fillna(0)

    # Anchor pEC50 mapping
    anchor_pec50 = {}
    for _, row in test_df.drop_duplicates("anchor_ligand_id").iterrows():
        idx = int(row["anchor_ligand_id"][1:])
        anchor_pec50[train_df.iloc[idx]["OCNT_ID"]] = row["best_train_pEC50"]

    # Phase 1 ground truth
    phase1_truth = dict(zip(phase1_df["Molecule Name"], phase1_df["pEC50"]))

    thresholds = [None, 3.0, 2.5, 2.0, 1.5, 1.0]
    rows = []

    for thresh in thresholds:
        preds = propagate(edge_df, anchor_pec50, max_error=thresh)
        # Evaluate against phase 1
        common = [(preds[c], phase1_truth[c])
                  for c in preds if c in phase1_truth]
        n_total = len(preds)
        n_phase1 = len(common)
        if n_phase1 >= 10:
            pred_vals = [x[0] for x in common]
            true_vals = [x[1] for x in common]
            rho, _ = spearmanr(true_vals, pred_vals)
            r, _ = pearsonr(true_vals, pred_vals)
            mae = np.mean(np.abs(np.array(true_vals) - np.array(pred_vals)))
        else:
            rho = r = mae = np.nan
        rows.append({
            "threshold": thresh if thresh is not None else "none",
            "compounds_connected": n_total,
            "phase1_evaluable": n_phase1,
            "spearman": rho,
            "pearson": r,
            "mae": mae,
        })
        t_str = f"{thresh}" if thresh is not None else "none"
        print(f"thresh={t_str:>5}  connected={n_total:3d}  "
              f"phase1_n={n_phase1:3d}  "
              f"rho={rho:.3f}  r={r:.3f}  MAE={mae:.3f}")

    sweep_df = pd.DataFrame(rows)
    out_path = outdir / "threshold_sweep.csv"
    sweep_df.to_csv(out_path, index=False)
    print(f"\nWrote {out_path}")

    # Recommend the threshold with best Spearman that keeps >=80% coverage
    valid = sweep_df[sweep_df["spearman"].notna()].copy()
    if len(valid):
        max_conn = valid["compounds_connected"].max()
        keep = valid[valid["compounds_connected"] >= 0.8 * max_conn]
        best = keep.loc[keep["spearman"].idxmax()]
        print(f"\nSuggested threshold (best Spearman with >=80% coverage):")
        print(f"  threshold={best['threshold']}  "
              f"connected={best['compounds_connected']}  "
              f"rho={best['spearman']:.3f}  MAE={best['mae']:.3f}")


if __name__ == "__main__":
    main()
