#!/usr/bin/env python3
"""
openfe/scripts/13_phase1_anchor_analysis.py

Analyze whether using unblinded Phase 1 test compounds as ADDITIONAL
anchors improves connectivity for Phase 2 (still-blind) compounds.

Phase 1 compounds (Analog Set 1) now have known pEC50. If they appear
as nodes in the RBFE networks (which they do, as OADMET test compounds),
they can serve as anchor points just like training compounds - giving
Phase 2 compounds more/closer reference points.

This script compares connectivity for Phase 2 compounds under two
anchor scenarios:
  A) Training anchors only (current approach)
  B) Training anchors + Phase 1 compounds with RBFE results

Outputs the gain: how many additional Phase 2 compounds become
connected, and how path lengths change.

Usage:
    python openfe/scripts/13_phase1_anchor_analysis.py \
        --edge-results openfe/all_edge_results.csv \
        --test-full openfe/test_full_with_clusters_and_anchors.csv \
        --train data/pxr-challenge_TRAIN.csv \
        --phase1 data/pxr-challenge_TEST_PHASE_1_UNBLINDED.csv \
        --outdir openfe
"""

import argparse
import re
from collections import defaultdict, deque
from pathlib import Path

import pandas as pd


def parse_edge_ligands(edge_name):
    body = edge_name[len("rbfe_"):] if edge_name.startswith("rbfe_") else edge_name
    parts = body.split("_")
    if len(parts) == 2:
        return parts[0], parts[1]
    return None, None


def build_graph(cluster_edges):
    graph = defaultdict(set)
    for _, row in cluster_edges.iterrows():
        if not row["both_done"]:
            continue
        a, b = parse_edge_ligands(row["edge"])
        if a is None:
            continue
        graph[a].add(b)
        graph[b].add(a)
    return graph


def reachable_from(graph, sources):
    """BFS from a set of source nodes; return {node: min_hops}."""
    dist = {s: 0 for s in sources if s in graph}
    queue = deque(dist.keys())
    while queue:
        node = queue.popleft()
        for nbr in graph[node]:
            if nbr not in dist:
                dist[nbr] = dist[node] + 1
                queue.append(nbr)
    return dist


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--edge-results", default="openfe/all_edge_results.csv")
    ap.add_argument("--test-full",
                    default="openfe/test_full_with_clusters_and_anchors.csv")
    ap.add_argument("--train", default="data/pxr-challenge_TRAIN.csv")
    ap.add_argument("--phase1",
                    default="data/pxr-challenge_TEST_PHASE_1_UNBLINDED.csv")
    ap.add_argument("--outdir", default="openfe")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    edge_df = pd.read_csv(args.edge_results)
    test_df = pd.read_csv(args.test_full)
    train_df = pd.read_csv(args.train)
    phase1_df = pd.read_csv(args.phase1)

    # Training anchors: OCNT_IDs with known pEC50
    train_anchor_ocnt = set()
    for _, row in test_df.drop_duplicates("anchor_ligand_id").iterrows():
        idx = int(row["anchor_ligand_id"][1:])
        train_anchor_ocnt.add(train_df.iloc[idx]["OCNT_ID"])

    # Phase 1 compounds (Molecule Name, these are OADMET names = graph nodes)
    phase1_names = set(phase1_df["Molecule Name"])

    # Map each compound to its set membership
    test_df = test_df.copy()
    test_df["set"] = test_df["Molecule Name"].apply(
        lambda x: 1 if x in phase1_names else 2)
    phase2_names = set(test_df[test_df["set"] == 2]["Molecule Name"])

    rows = []
    summary = {
        "phase2_connected_train_only": 0,
        "phase2_connected_with_phase1": 0,
        "phase2_total_in_networks": 0,
    }

    for cluster_id, cluster_edges in edge_df.groupby("cluster_id"):
        graph = build_graph(cluster_edges)
        if not graph:
            continue

        nodes = set(graph.keys())
        train_anchors = {n for n in nodes if n in train_anchor_ocnt}
        phase1_in_cluster = {n for n in nodes if n in phase1_names}
        phase2_in_cluster = {n for n in nodes if n in phase2_names}

        if not phase2_in_cluster:
            continue
        summary["phase2_total_in_networks"] += len(phase2_in_cluster)

        # Scenario A: training anchors only
        reach_A = reachable_from(graph, train_anchors)
        # Scenario B: training anchors + phase 1 compounds
        reach_B = reachable_from(graph, train_anchors | phase1_in_cluster)

        for comp in phase2_in_cluster:
            conn_A = comp in reach_A
            conn_B = comp in reach_B
            hops_A = reach_A.get(comp)
            hops_B = reach_B.get(comp)
            if conn_A:
                summary["phase2_connected_train_only"] += 1
            if conn_B:
                summary["phase2_connected_with_phase1"] += 1
            rows.append({
                "cluster_id": cluster_id,
                "compound": comp,
                "connected_train_only": conn_A,
                "connected_with_phase1": conn_B,
                "hops_train_only": hops_A,
                "hops_with_phase1": hops_B,
                "newly_connected": conn_B and not conn_A,
                "shorter_path": (conn_A and conn_B and hops_B < hops_A),
                "n_train_anchors_in_cluster": len(train_anchors),
                "n_phase1_in_cluster": len(phase1_in_cluster),
            })

    result_df = pd.DataFrame(rows)
    out_path = outdir / "phase1_anchor_analysis.csv"
    result_df.to_csv(out_path, index=False)

    print("=" * 60)
    print("PHASE 1 AS ADDITIONAL ANCHORS - IMPACT ON PHASE 2")
    print("=" * 60)
    print(f"Phase 2 compounds in RBFE networks: "
          f"{summary['phase2_total_in_networks']}")
    print(f"\nConnected with TRAINING anchors only: "
          f"{summary['phase2_connected_train_only']}")
    print(f"Connected with TRAINING + PHASE 1 anchors: "
          f"{summary['phase2_connected_with_phase1']}")
    gain = (summary["phase2_connected_with_phase1"]
            - summary["phase2_connected_train_only"])
    print(f"\nNewly connected Phase 2 compounds: {gain}")

    if len(result_df):
        newly = result_df[result_df["newly_connected"]]
        shorter = result_df[result_df["shorter_path"]]
        print(f"Phase 2 compounds with NEW connection: {len(newly)}")
        print(f"Phase 2 compounds with SHORTER path: {len(shorter)}")

        if len(newly):
            print(f"\nNewly connected compounds by cluster:")
            print(newly.groupby("cluster_id").size().to_string())

    print(f"\nWrote {out_path}")


if __name__ == "__main__":
    main()
