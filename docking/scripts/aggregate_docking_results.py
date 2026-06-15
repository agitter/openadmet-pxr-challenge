#!/usr/bin/env python3
"""
docking/scripts/aggregate_docking_results.py

Aggregate GNINA ensemble-docking results (124 cluster representatives x
~62 receptors, 5 poses/receptor) into analysis-ready tables. Intended
to be run LOCALLY (uses RDKit).

Reads:
  results/<cluster_id>_summary.csv
      one row per (cluster_id, pdb_id): status in
      {success, FAILED, FAILED_NO_RECEPTOR}
  results/<cluster_id>/<cluster_id>__<pdb_id>_docked.sdf.gz
      gnina output, up to --num-modes poses per (cluster, receptor),
      each tagged with minimizedAffinity, CNNscore, CNNaffinity.

Writes (to --outdir):
  docking_poses_long.csv
      One row per (cluster_id, pdb_id, pose_rank):
      minimizedAffinity, CNNscore, CNNaffinity.
      pose_rank is the 1-based order of appearance in the SDF (gnina's
      internal ordering - NOT necessarily sorted by CNNscore).

  docking_receptor_best.csv
      One row per (cluster_id, pdb_id): the pose with the highest
      CNNscore for that receptor, i.e. gnina's most confident binding
      mode in that pocket conformation.

  docking_cluster_summary.csv
      One row per cluster_id:
        - best (receptor, pose) overall by CNNscore
        - mean/std of best-per-receptor CNNscore across all receptors
          (a consistency signal: low std = most pocket conformations
          agree this is a confident pose; high std = only a few
          receptors support a confident binding mode)
        - top-N receptors by CNNscore (candidates for RBFE starting
          structures / ensemble templates)
        - receptor success/failure counts from the summary CSVs

Usage:
    pip install rdkit pandas
    python aggregate_docking_results.py \
        --results-dir results \
        --outdir docking_analysis \
        --top-n 3
"""

import argparse
import gzip
from pathlib import Path

import numpy as np
import pandas as pd

from rdkit import Chem
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")


def parse_sdf_gz(path):
    """Parse a gzipped multi-record SDF via RDKit, returning a list of
    SD-tag dicts (one per pose, in file order -> pose_rank). A failed
    molblock (mol is None) yields an empty dict so rank ordering of
    the remaining poses is preserved."""
    poses = []
    with gzip.open(path, "rb") as f:
        supplier = Chem.ForwardSDMolSupplier(f, sanitize=False, removeHs=False)
        for mol in supplier:
            if mol is None:
                poses.append({})
                continue
            poses.append(mol.GetPropsAsDict(includePrivate=False,
                                             includeComputed=False))
    return poses


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results-dir", default="results")
    ap.add_argument("--outdir", default="docking_analysis")
    ap.add_argument("--top-n", type=int, default=3,
                    help="Number of top receptors (by CNNscore) to "
                         "record per cluster (default: 3)")
    ap.add_argument("--id-pattern", default=None,
                    help="If set, only include results whose cluster_id "
                         "(as a string) matches this regex (anchored at "
                         "the start, via str.match). Useful when "
                         "results/ contains multiple campaigns with "
                         "different ID schemes - e.g. '^[TA]\\d+$' to "
                         "select only the extended T###/A#### campaign "
                         "and exclude the original integer cluster_id "
                         "directories. Default: no filtering (use all).")
    args = ap.parse_args()

    results_dir = Path(args.results_dir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # -----------------------------------------------------------
    # 1. Read all per-cluster summary CSVs -> overall status table
    # -----------------------------------------------------------
    summary_files = sorted(results_dir.glob("*_summary.csv"))
    print(f"Found {len(summary_files)} summary CSVs")
    if not summary_files:
        raise SystemExit(f"No *_summary.csv files found in {results_dir}")

    status_df = pd.concat(
        (pd.read_csv(f) for f in summary_files), ignore_index=True)
    print(f"Total (cluster, receptor) pairs (before filtering): {len(status_df)}")

    if args.id_pattern:
        mask = status_df["cluster_id"].astype(str).str.match(args.id_pattern)
        n_dropped = (~mask).sum()
        status_df = status_df[mask].reset_index(drop=True)
        print(f"Filtered to cluster_id matching '{args.id_pattern}': "
              f"{len(status_df)} pairs kept, {n_dropped} dropped")

    print(f"Total (cluster, receptor) pairs: {len(status_df)}")
    print(status_df["status"].value_counts().to_string())

    n_clusters = status_df["cluster_id"].nunique()
    print(f"Clusters represented: {n_clusters}")

    # -----------------------------------------------------------
    # 2. Parse SDFs for every successful (cluster, receptor) pair
    # -----------------------------------------------------------
    long_rows = []
    n_parsed, n_missing = 0, 0
    for _, row in status_df[status_df["status"] == "success"].iterrows():
        cluster_id, pdb_id, ligand_name = (
            row["cluster_id"], row["pdb_id"], row["ligand_name"])
        sdf_path = (results_dir / str(cluster_id)
                     / f"{cluster_id}__{pdb_id}_docked.sdf.gz")
        if not sdf_path.exists():
            n_missing += 1
            continue
        for rank, props in enumerate(parse_sdf_gz(sdf_path), start=1):
            long_rows.append({
                "cluster_id": cluster_id,
                "ligand_name": ligand_name,
                "pdb_id": pdb_id,
                "pose_rank": rank,
                "minimizedAffinity": props.get("minimizedAffinity", np.nan),
                "CNNscore": props.get("CNNscore", np.nan),
                "CNNaffinity": props.get("CNNaffinity", np.nan),
            })
        n_parsed += 1

    print(f"\nParsed {n_parsed} SDFs "
          f"({n_missing} expected-success SDFs missing on disk)")

    long_df = pd.DataFrame(long_rows)
    for col in ("minimizedAffinity", "CNNscore", "CNNaffinity"):
        long_df[col] = pd.to_numeric(long_df[col], errors="coerce")

    long_df.to_csv(outdir / "docking_poses_long.csv", index=False)
    print(f"Wrote {outdir/'docking_poses_long.csv'} "
          f"({len(long_df)} pose rows)")

    if long_df.empty:
        print("No poses parsed - stopping before aggregation steps.")
        return

    # -----------------------------------------------------------
    # 3. Per (cluster, receptor): pose with the highest CNNscore
    # -----------------------------------------------------------
    idx = long_df.groupby(["cluster_id", "pdb_id"])["CNNscore"].idxmax()
    receptor_best = (long_df.loc[idx]
                      .sort_values(["cluster_id", "CNNscore"],
                                   ascending=[True, False])
                      .reset_index(drop=True))
    receptor_best.to_csv(outdir / "docking_receptor_best.csv", index=False)
    print(f"Wrote {outdir/'docking_receptor_best.csv'} "
          f"({len(receptor_best)} cluster x receptor rows)")

    # -----------------------------------------------------------
    # 4. Per cluster: overall best (receptor, pose) + consistency
    # -----------------------------------------------------------
    cluster_rows = []
    for cluster_id, group in receptor_best.groupby("cluster_id"):
        best = group.loc[group["CNNscore"].idxmax()]

        cluster_status = status_df[status_df["cluster_id"] == cluster_id]
        n_success = (cluster_status["status"] == "success").sum()
        n_total = len(cluster_status)

        top_n = (group.sort_values("CNNscore", ascending=False)
                 .head(args.top_n))

        cluster_rows.append({
            "cluster_id": cluster_id,
            "ligand_name": best["ligand_name"],
            "best_pdb_id": best["pdb_id"],
            "best_pose_rank": int(best["pose_rank"]),
            "best_CNNscore": best["CNNscore"],
            "best_CNNaffinity": best["CNNaffinity"],
            "best_minimizedAffinity": best["minimizedAffinity"],
            "mean_CNNscore_across_receptors": group["CNNscore"].mean(),
            "std_CNNscore_across_receptors": group["CNNscore"].std(),
            f"top{args.top_n}_receptors": ",".join(top_n["pdb_id"].tolist()),
            f"top{args.top_n}_CNNscores": ",".join(
                f"{v:.3f}" for v in top_n["CNNscore"].tolist()),
            "n_receptors_success": n_success,
            "n_receptors_total": n_total,
        })

    cluster_summary = (pd.DataFrame(cluster_rows)
                        .sort_values("cluster_id")
                        .reset_index(drop=True))
    cluster_summary.to_csv(outdir / "docking_cluster_summary.csv", index=False)
    print(f"Wrote {outdir/'docking_cluster_summary.csv'} "
          f"({len(cluster_summary)} clusters)")

    # -----------------------------------------------------------
    # 5. Quick report
    # -----------------------------------------------------------
    print("\n=== Quick stats ===")
    print(f"Clusters with results: {len(cluster_summary)} / 124 expected")
    print(f"Best CNNscore distribution:")
    print(cluster_summary["best_CNNscore"].describe().to_string())

    complete = (cluster_summary["n_receptors_success"]
                 == cluster_summary["n_receptors_total"])
    print(f"\nClusters with full receptor coverage: "
          f"{complete.sum()} / {len(cluster_summary)}")

    incomplete = cluster_summary[~complete]
    if len(incomplete):
        print("\nClusters with incomplete receptor coverage:")
        print(incomplete[["cluster_id", "n_receptors_success",
                           "n_receptors_total"]].to_string(index=False))

    low_conf = cluster_summary[cluster_summary["best_CNNscore"] < 0.5]
    if len(low_conf):
        print(f"\n{len(low_conf)} clusters with best_CNNscore < 0.5 "
              f"(low pose confidence across ALL receptors - may need "
              f"closer inspection / Boltz-2 cross-check):")
        print(low_conf[["cluster_id", "ligand_name", "best_pdb_id",
                         "best_CNNscore"]].to_string(index=False))


if __name__ == "__main__":
    main()
