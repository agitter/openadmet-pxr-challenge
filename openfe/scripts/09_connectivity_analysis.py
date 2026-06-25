#!/usr/bin/env python3
"""
openfe/scripts/09_connectivity_analysis.py

Determine the real impact of failed edges on compound coverage.

For each cluster, build a graph using ONLY the successfully completed
edges (both legs done). Then determine which test compounds are still
connected to at least one training anchor (and thus have an RBFE-based
pEC50 prediction via the edge path) versus which are disconnected
(and would need a docking fallback or edge salvage).

This tells us whether salvaging failed edges is worth it:
  - If few compounds are disconnected, failed edges are mostly redundant
  - If many are disconnected, salvaging specific edges has high value

Inputs:
  - openfe/final_edge_report.csv (from 08_final_report.py)
  - openfe/test_full_with_clusters_and_anchors.csv (compound->cluster->anchor)

Outputs:
  - openfe/connectivity_report.csv (per test compound: connected?, n_paths)
  - openfe/salvage_priority_edges.csv (failed edges ranked by how many
    compounds they would reconnect)

Usage:
    python openfe/scripts/09_connectivity_analysis.py \
        --edge-report openfe/final_edge_report.csv \
        --test-full openfe/test_full_with_clusters_and_anchors.csv \
        --outdir openfe
"""

import argparse
from collections import defaultdict
from pathlib import Path

import pandas as pd


def parse_edge_ligands(edge_name):
    """rbfe_<ligA>_<ligB> -> (ligA, ligB). Ligand names contain
    hyphens (OADMET-0006503, OCNT-2317296) but no underscores, so
    split on underscore after stripping the rbfe_ prefix."""
    body = edge_name[len("rbfe_"):] if edge_name.startswith("rbfe_") else edge_name
    parts = body.split("_")
    # Ligand names have no underscores, so exactly 2 parts
    if len(parts) == 2:
        return parts[0], parts[1]
    return None, None


def build_graph(edges):
    """Adjacency list from a list of (ligA, ligB) tuples."""
    graph = defaultdict(set)
    for a, b in edges:
        graph[a].add(b)
        graph[b].add(a)
    return graph


def connected_component(graph, start):
    """BFS from start, return set of reachable nodes."""
    seen = {start}
    queue = [start]
    while queue:
        node = queue.pop()
        for nbr in graph[node]:
            if nbr not in seen:
                seen.add(nbr)
                queue.append(nbr)
    return seen


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--edge-report", default="openfe/final_edge_report.csv")
    ap.add_argument("--test-full",
                    default="openfe/test_full_with_clusters_and_anchors.csv")
    ap.add_argument("--outdir", default="openfe")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    edge_df = pd.read_csv(args.edge_report)
    test_full = pd.read_csv(args.test_full)

    # Map: cluster_id -> set of anchor ligand_ids (the OCNT names),
    # and cluster_id -> set of test compound names (OADMET names)
    # The edge report uses Molecule Name / OCNT_ID which appear in edges.
    # test_full has: Molecule Name (test), anchor info.

    # Build per-cluster anchor name set from test_full.
    # Anchors are training compounds; their name in edges is the OCNT_ID.
    # We need a mapping from anchor_ligand_id (A####) to OCNT_ID.
    # test_full has best_train_OCNT or similar - check columns.
    # We'll treat any node starting with OCNT as an anchor.

    conn_rows = []
    salvage_value = defaultdict(int)  # (cluster, edge) -> compounds reconnected

    for cluster_id, grp in edge_df.groupby("cluster_id"):
        # Working edges: both legs done
        working = []
        failed = []
        for _, row in grp.iterrows():
            a, b = parse_edge_ligands(row["edge"])
            if a is None:
                continue
            if row["both_done"]:
                working.append((a, b))
            else:
                failed.append((a, b, row["edge"]))

        graph = build_graph(working)

        # Anchors in this cluster = nodes starting with OCNT
        all_nodes = set()
        for a, b in working:
            all_nodes.add(a)
            all_nodes.add(b)
        for a, b, _ in failed:
            all_nodes.add(a)
            all_nodes.add(b)
        anchors = {n for n in all_nodes if n.startswith("OCNT")}
        test_compounds = {n for n in all_nodes if n.startswith("OADMET")}

        # Which test compounds are connected to an anchor via working edges?
        anchor_reachable = set()
        for anchor in anchors:
            if anchor in graph:
                anchor_reachable |= connected_component(graph, anchor)

        for tc in test_compounds:
            connected = tc in anchor_reachable
            conn_rows.append({
                "cluster_id": cluster_id,
                "compound": tc,
                "connected_to_anchor": connected,
            })

        # Salvage value: for each failed edge, how many disconnected
        # test compounds would it reconnect if it worked?
        disconnected = test_compounds - anchor_reachable
        if disconnected and failed:
            # Try adding each failed edge and see what reconnects
            for a, b, ename in failed:
                test_graph = build_graph(working + [(a, b)])
                reachable = set()
                for anchor in anchors:
                    if anchor in test_graph:
                        reachable |= connected_component(test_graph, anchor)
                newly = (disconnected & reachable)
                salvage_value[(cluster_id, ename)] = len(newly)

    conn_df = pd.DataFrame(conn_rows)
    conn_path = outdir / "connectivity_report.csv"
    conn_df.to_csv(conn_path, index=False)

    n_connected = conn_df["connected_to_anchor"].sum()
    n_total = len(conn_df)
    print("=" * 60)
    print("CONNECTIVITY ANALYSIS")
    print("=" * 60)
    print(f"Test compounds in RBFE networks: {n_total}")
    print(f"Connected to anchor (RBFE-predictable): {n_connected} "
          f"({100*n_connected/n_total:.1f}%)")
    print(f"Disconnected (need salvage or docking): {n_total - n_connected} "
          f"({100*(n_total-n_connected)/n_total:.1f}%)")

    # Salvage priority
    salvage_rows = [
        {"cluster_id": k[0], "edge": k[1], "compounds_reconnected": v}
        for k, v in salvage_value.items() if v > 0
    ]
    if salvage_rows:
        salvage_df = pd.DataFrame(salvage_rows).sort_values(
            "compounds_reconnected", ascending=False)
        salvage_path = outdir / "salvage_priority_edges.csv"
        salvage_df.to_csv(salvage_path, index=False)
        print(f"\nSALVAGE PRIORITY:")
        print(f"Failed edges that would reconnect >=1 compound: "
              f"{len(salvage_df)}")
        print(f"Total reconnectable compounds (if all salvaged): "
              f"{salvage_df['compounds_reconnected'].sum()}")
        print(f"\nTop 15 highest-value edges to salvage:")
        print(salvage_df.head(15).to_string(index=False))
        print(f"\nWrote {salvage_path}")
    else:
        print("\nNo single-edge salvage reconnects compounds "
              "(failures may be redundant or need multiple edges).")

    print(f"\nWrote {conn_path}")


if __name__ == "__main__":
    main()
