#!/usr/bin/env python3
"""
openfe/scripts/12_propagate_pec50.py

Propagate RBFE ddG values along MST networks to predict absolute pEC50
for test compounds connected to an anchor.

This script does ONE thing: build per-cluster graphs from completed
edges, compute the ddG path from anchor to each connected test
compound, and convert to a raw (uncalibrated) pEC50.

    pEC50(test) = pEC50(anchor) - ddG_path / RT_ln10

where RT_ln10 = 1.3638 kcal/mol at 298K.

Sign convention:
    ddG = dG_complex - dG_solvent  (per edge, directed A->B)
    More negative ddG_path = stronger binding = higher pEC50.

Calibration, docking fallback, and final submission formatting are
handled by SEPARATE downstream scripts.

Output:
    openfe/rbfe_predictions.csv - one row per connected test compound:
        Molecule Name, cluster_id, anchor_ocnt, anchor_pEC50,
        path_ddg_kcal_mol, n_hops, path, pred_pEC50_raw

Usage:
    python openfe/scripts/12_propagate_pec50.py \
        --edge-results openfe/all_edge_results.csv \
        --test-full openfe/test_full_with_clusters_and_anchors.csv \
        --train data/pxr-challenge_TRAIN.csv \
        --outdir openfe
"""

import argparse
import math
import re
from collections import defaultdict, deque
from pathlib import Path

import pandas as pd

# RT * ln(10) at 298K in kcal/mol: 0.0019872 * 298.15 * ln(10)
RT_LN10 = 1.3638


def parse_edge_ligands(edge_name):
    """rbfe_<ligA>_<ligB> -> (ligA, ligB). Ligand names contain hyphens
    but no underscores, so split on underscore after the rbfe_ prefix."""
    body = edge_name[len("rbfe_"):] if edge_name.startswith("rbfe_") else edge_name
    parts = body.split("_")
    if len(parts) == 2:
        return parts[0], parts[1]
    return None, None


def build_graph(cluster_edges, edge_error=None, edge_overlap=None):
    """Directed-ddG adjacency: node -> {neighbor: (ddG, edge_err, overlap)}.
    A->B stores ddG; B->A stores -ddG so a path sum is consistent.
    edge_err is the propagated ddG MBAR error for that edge; overlap is
    the min MBAR overlap scalar across the edge's two legs. Both are the
    same in either direction."""
    graph = defaultdict(dict)
    for _, row in cluster_edges.iterrows():
        if not row["both_done"]:
            continue
        a, b = parse_edge_ligands(row["edge"])
        if a is None:
            continue
        ddg = row["ddg"]
        err = edge_error.get(row["edge"]) if edge_error is not None else None
        ovl = edge_overlap.get(row["edge"]) if edge_overlap is not None else None
        graph[a][b] = (ddg, err, ovl)
        graph[b][a] = (-ddg, err, ovl)
    return graph


def bfs_paths(graph, start):
    """BFS from start; fewest-hops path to each node (unique in an MST).
    Returns {node: (cumulative_ddg, n_hops, path, path_err_quad,
    max_edge_err, min_overlap)} where path_err_quad is the sqrt-sum-of-
    squares of per-edge errors, max_edge_err is the worst single edge
    error, and min_overlap is the smallest MBAR overlap on the path."""
    visited = {start: (0.0, 0, [start], 0.0, 0.0, float("inf"))}
    queue = deque([start])
    while queue:
        node = queue.popleft()
        (curr_ddg, curr_hops, curr_path, curr_var,
         curr_max, curr_minov) = visited[node]
        for neighbor, (ddg, err, ovl) in graph[node].items():
            if neighbor not in visited:
                e = err if (err is not None and not math.isnan(err)) else 0.0
                new_var = curr_var + e * e
                new_max = max(curr_max, e)
                if ovl is not None and not math.isnan(ovl):
                    new_minov = min(curr_minov, ovl)
                else:
                    new_minov = curr_minov
                visited[neighbor] = (
                    curr_ddg + ddg, curr_hops + 1,
                    curr_path + [neighbor], new_var, new_max, new_minov)
                queue.append(neighbor)
    out = {}
    for n, (d, h, p, v, m, mo) in visited.items():
        min_ov = mo if mo != float("inf") else None
        out[n] = (d, h, p, math.sqrt(v), m, min_ov)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--edge-results", default="openfe/all_edge_results.csv")
    ap.add_argument("--edge-convergence",
                    default="openfe/edge_convergence.csv",
                    help="Per-edge MBAR errors (from 15_convergence_analysis). "
                         "Used to compute path_ddg_error and "
                         "max_edge_error_on_path. Optional.")
    ap.add_argument("--test-full",
                    default="openfe/test_full_with_clusters_and_anchors.csv")
    ap.add_argument("--train", default="data/pxr-challenge_TRAIN.csv")
    ap.add_argument("--outdir", default="openfe")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    edge_df = pd.read_csv(args.edge_results)
    test_df = pd.read_csv(args.test_full)
    train_df = pd.read_csv(args.train)

    # Build per-edge ddG error lookup: sqrt(complex_err^2 + solvent_err^2)
    # and per-edge min overlap (min of the two legs' overlap scalars).
    edge_error = {}
    edge_overlap = {}
    conv_path = Path(args.edge_convergence)
    if conv_path.exists():
        conv = pd.read_csv(conv_path)
        ecols = [c for c in conv.columns if c.endswith("_dG_error")]
        ocols = [c for c in conv.columns if c.endswith("_min_overlap")]
        for _, r in conv.iterrows():
            errs = [r[c] for c in ecols if pd.notna(r[c])]
            if errs:
                edge_error[r["edge"]] = math.sqrt(sum(e * e for e in errs))
            ovls = [r[c] for c in ocols if pd.notna(r[c])]
            if ovls:
                edge_overlap[r["edge"]] = min(ovls)
        print(f"Loaded per-edge errors for {len(edge_error)} edges, "
              f"overlaps for {len(edge_overlap)} edges")
    else:
        print(f"No edge-convergence file at {conv_path}; "
              f"path error/overlap columns will be 0/None.")

    # Build anchor OCNT_ID -> pEC50 mapping.
    # anchor_ligand_id is "A####" = row index into the training set.
    anchor_pec50 = {}
    for _, row in test_df.drop_duplicates("anchor_ligand_id").iterrows():
        idx = int(row["anchor_ligand_id"][1:])
        train_row = train_df.iloc[idx]
        anchor_pec50[train_row["OCNT_ID"]] = row["best_train_pEC50"]
    print(f"Anchor pEC50 mappings: {len(anchor_pec50)}")

    rows = []
    n_clusters_with_anchor = 0

    for cluster_id, cluster_edges in edge_df.groupby("cluster_id"):
        graph = build_graph(cluster_edges, edge_error=edge_error,
                            edge_overlap=edge_overlap)
        if not graph:
            continue

        all_nodes = set(graph.keys())
        anchors = {n for n in all_nodes
                   if n.startswith("OCNT") and n in anchor_pec50}
        test_compounds = {n for n in all_nodes if n.startswith("OADMET")}
        if not anchors:
            continue
        n_clusters_with_anchor += 1

        # For each test compound, take the nearest anchor (fewest hops)
        best = {}  # compound -> (path_ddg, n_hops, anchor, path, perr, pmax, pmo)
        for anchor in anchors:
            paths = bfs_paths(graph, anchor)
            for compound in test_compounds:
                if compound in paths:
                    path_ddg, n_hops, path, perr, pmax, pmo = paths[compound]
                    if compound not in best or n_hops < best[compound][1]:
                        best[compound] = (path_ddg, n_hops, anchor, path,
                                          perr, pmax, pmo)

        for compound, (path_ddg, n_hops, anchor, path, perr, pmax, pmo) in \
                best.items():
            pec50_anchor = anchor_pec50[anchor]
            pred = pec50_anchor - path_ddg / RT_LN10
            rows.append({
                "Molecule Name": compound,
                "cluster_id": cluster_id,
                "anchor_ocnt": anchor,
                "anchor_pEC50": pec50_anchor,
                "path_ddg_kcal_mol": path_ddg,
                "n_hops": n_hops,
                "path": " -> ".join(path),
                "pred_pEC50_raw": pred,
                "path_ddg_error": perr,
                "max_edge_error_on_path": pmax,
                "min_overlap_on_path": pmo,
            })

    rbfe_df = pd.DataFrame(rows)
    out_path = outdir / "rbfe_predictions.csv"
    rbfe_df.to_csv(out_path, index=False)

    print(f"Clusters with at least one anchor: {n_clusters_with_anchor}")
    print(f"RBFE-connected test compounds: {len(rbfe_df)}")
    if len(rbfe_df):
        print(f"\nRaw pEC50 prediction distribution:")
        print(f"  mean={rbfe_df['pred_pEC50_raw'].mean():.2f}  "
              f"median={rbfe_df['pred_pEC50_raw'].median():.2f}  "
              f"min={rbfe_df['pred_pEC50_raw'].min():.2f}  "
              f"max={rbfe_df['pred_pEC50_raw'].max():.2f}")
        print(f"\nPath length distribution (hops from anchor):")
        print(rbfe_df["n_hops"].value_counts().sort_index().to_string())
        print(f"\nddG path magnitude:")
        print(f"  mean |ddG_path|={rbfe_df['path_ddg_kcal_mol'].abs().mean():.2f}  "
              f"max |ddG_path|={rbfe_df['path_ddg_kcal_mol'].abs().max():.2f}")
    print(f"\nWrote {out_path}")


if __name__ == "__main__":
    main()
