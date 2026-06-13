#!/usr/bin/env python3
"""
analysis/01_cluster_test_set.py

ECFP4/Tanimoto clustering of the 513-compound PXR challenge test set,
matching the dataset's own construction criterion (analogs selected at
ECFP4 Tanimoto > 0.4 to 63 confirmed hits).

Produces:
  - test_with_clusters.csv : the test set with an added 'cluster_id'
    column (single-linkage clusters at the chosen Tanimoto threshold)
  - cluster_summary.csv : one row per cluster (size, best-similarity-
    based representative compound)
  - test_tanimoto_sim.npy : full 513x513 pairwise Tanimoto similarity
    matrix (ECFP4, 2048 bits, radius 2) - reusable for other thresholds
    or downstream analyses without recomputation

Usage:
    pip install rdkit pandas numpy
    python 01_cluster_test_set.py \
        --test pxr-challenge_TEST_BLINDED.csv \
        --threshold 0.5
"""

import argparse
import pickle
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")


def get_fp(smi, radius=2, nbits=2048):
    if pd.isna(smi):
        return None
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    return AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=nbits)


def cluster_at_threshold(fps, threshold):
    """Single-linkage clustering via union-find on a Tanimoto similarity
    matrix. Returns (sim_matrix, cluster_id_per_compound)."""
    n = len(fps)
    sim = np.zeros((n, n))
    for i in range(n):
        sim[i, :] = DataStructs.BulkTanimotoSimilarity(fps[i], fps)

    adj = sim >= threshold
    np.fill_diagonal(adj, False)

    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            x = parent[x]
        return x

    def union(x, y):
        px, py = find(x), find(y)
        if px != py:
            parent[px] = py

    for i in range(n):
        for j in range(i + 1, n):
            if adj[i, j]:
                union(i, j)

    cluster_id = [find(i) for i in range(n)]
    return sim, cluster_id


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--test", required=True,
                    help="Path to pxr-challenge_TEST_BLINDED.csv")
    ap.add_argument("--threshold", type=float, default=0.5,
                    help="Tanimoto similarity threshold for clustering "
                         "(default: 0.5). The dataset itself was built "
                         "using >0.4 similarity to 63 hits; 0.5 gives "
                         "tighter, more RBFE-suitable clusters.")
    ap.add_argument("--outdir", default=".", help="Output directory")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    test = pd.read_csv(args.test)
    print(f"Loaded {len(test)} test compounds")

    test["fp"] = test["SMILES"].apply(get_fp)
    n_valid = test["fp"].notna().sum()
    print(f"Valid fingerprints: {n_valid} / {len(test)}")
    if n_valid != len(test):
        print("WARNING: some SMILES failed to parse - check input data")
        test = test[test["fp"].notna()].reset_index(drop=True)

    fps = test["fp"].tolist()
    sim, cluster_id = cluster_at_threshold(fps, args.threshold)
    test["cluster_id"] = cluster_id

    n_clusters = len(set(cluster_id))
    sizes = test["cluster_id"].value_counts()
    print(f"\nThreshold {args.threshold}: {n_clusters} clusters")
    print(f"Singletons: {(sizes == 1).sum()}")
    print(f"Largest 10 cluster sizes: {sizes.head(10).tolist()}")
    print(f"Compounds in clusters of size>=2: {sizes[sizes >= 2].sum()}")

    # Build cluster summary: pick the compound closest to the cluster
    # centroid (highest mean similarity to other cluster members) as
    # the representative.
    summary_rows = []
    for cid, group in test.groupby("cluster_id"):
        idxs = group.index.tolist()
        if len(idxs) == 1:
            rep_idx = idxs[0]
        else:
            # mean similarity to other members of the same cluster
            mean_sims = []
            for i in idxs:
                sims = sim[i, idxs]
                mean_sims.append(sims.sum() / (len(idxs) - 1) if len(idxs) > 1 else 0)
            rep_idx = idxs[int(np.argmax(mean_sims))]
        summary_rows.append({
            "cluster_id": cid,
            "cluster_size": len(idxs),
            "representative_test_compound": test.loc[rep_idx, "Molecule Name"],
            "representative_smiles": test.loc[rep_idx, "SMILES"],
        })

    cluster_summary = pd.DataFrame(summary_rows).sort_values(
        "cluster_size", ascending=False)

    # Save outputs
    test.drop(columns=["fp"]).to_csv(outdir / "test_with_clusters.csv", index=False)
    cluster_summary.to_csv(outdir / "cluster_summary.csv", index=False)
    np.save(outdir / "test_tanimoto_sim.npy", sim)
    with open(outdir / "test_fps.pkl", "wb") as f:
        pickle.dump({"molecule_names": test["Molecule Name"].tolist(),
                      "cluster_id": cluster_id,
                      "threshold": args.threshold}, f)

    print(f"\nWrote:")
    print(f"  {outdir / 'test_with_clusters.csv'}")
    print(f"  {outdir / 'cluster_summary.csv'}")
    print(f"  {outdir / 'test_tanimoto_sim.npy'}")
    print(f"  {outdir / 'test_fps.pkl'}")


if __name__ == "__main__":
    main()
