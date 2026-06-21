#!/usr/bin/env python3
"""
openfe/scripts/05_setup_production.py

Create per-job directories for the production openfe quickrun campaign.
Each transformation JSON from the network planning results gets its own
directory under production/<cluster_id>/<transform_name>/, containing:
  - The transformation JSON (copied)
  - network_setup.json (copied from the cluster's network planning output)
  - quickrun_output/ (empty directory, will hold checkpoints during run)

Also writes production/transform_list.txt as the HTCondor queue source.

Usage:
    python openfe/scripts/05_setup_production.py \
        --network-results openfe/results \
        --outdir openfe/production
"""

import argparse
import shutil
from pathlib import Path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--network-results", default="openfe/results",
                    help="Directory containing network_setup_<cluster_id>/ "
                         "subdirectories from planning jobs")
    ap.add_argument("--outdir", default="openfe/production")
    args = ap.parse_args()

    network_results = Path(args.network_results)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Find all cluster directories
    cluster_dirs = sorted(
        d for d in network_results.iterdir()
        if d.is_dir() and d.name.startswith("network_setup_")
    )
    print(f"Found {len(cluster_dirs)} cluster network directories")

    transform_entries = []
    n_total = 0
    n_skipped = 0

    for cluster_dir in cluster_dirs:
        cluster_id = cluster_dir.name.replace("network_setup_", "")
        network_json = cluster_dir / "network_setup.json"
        transforms_dir = cluster_dir / "transformations"

        if not network_json.exists():
            print(f"  WARNING: {cluster_id} has no network_setup.json, "
                  f"skipping")
            continue
        if not transforms_dir.exists():
            print(f"  WARNING: {cluster_id} has no transformations/, "
                  f"skipping")
            continue

        json_files = sorted(transforms_dir.glob("*.json"))
        for json_file in json_files:
            transform_name = json_file.stem  # filename without .json
            job_dir = outdir / cluster_id / transform_name
            job_dir.mkdir(parents=True, exist_ok=True)

            # Copy transformation JSON
            dest_json = job_dir / json_file.name
            if not dest_json.exists():
                shutil.copy(json_file, dest_json)

            # Copy network_setup.json
            dest_network = job_dir / "network_setup.json"
            if not dest_network.exists():
                shutil.copy(network_json, dest_network)

            # Create empty quickrun_output directory for checkpoints
            qo = job_dir / "quickrun_output"
            qo.mkdir(exist_ok=True)

            transform_entries.append(f"{cluster_id}, {transform_name}")
            n_total += 1

    # Write queue source file
    list_path = outdir / "transform_list.txt"
    with open(list_path, "w") as f:
        f.write("\n".join(transform_entries) + "\n")

    print(f"\nTotal transformation jobs: {n_total}")
    print(f"Skipped: {n_skipped}")
    print(f"Wrote {list_path}")
    print(f"\nJob directories created under {outdir}/")
    print(f"Each contains: <transform>.json, network_setup.json, "
          f"quickrun_output/")


if __name__ == "__main__":
    main()
