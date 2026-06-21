#!/usr/bin/env python3
"""
openfe/scripts/05_setup_production.py

Create per-job directories for the production openfe quickrun campaign.
Each transformation JSON gets its own directory under
production/<cluster_id>/<transform_name>/, used as the HTCondor
initialdir for that job.

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
    ap.add_argument("--script", default="openfe/scripts/run_quickrun.sh",
                    help="Path to run_quickrun.sh (copied into each "
                         "job directory for reliable HTCondor transfer)")
    ap.add_argument("--outdir", default="openfe/production")
    args = ap.parse_args()

    network_results = Path(args.network_results)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    cluster_dirs = sorted(
        d for d in network_results.iterdir()
        if d.is_dir() and d.name.startswith("network_setup_")
    )
    print(f"Found {len(cluster_dirs)} cluster network directories")

    transform_entries = []
    n_total = 0

    for cluster_dir in cluster_dirs:
        cluster_id = cluster_dir.name.replace("network_setup_", "")
        network_json = cluster_dir / "network_setup.json"
        transforms_dir = cluster_dir / "transformations"

        if not network_json.exists():
            print(f"  WARNING: {cluster_id} missing network_setup.json")
            continue
        if not transforms_dir.exists():
            print(f"  WARNING: {cluster_id} missing transformations/")
            continue

        for json_file in sorted(transforms_dir.glob("*.json")):
            transform_name = json_file.stem
            job_dir = outdir / cluster_id / transform_name
            job_dir.mkdir(parents=True, exist_ok=True)

            # Copy transformation JSON into job dir
            dest_json = job_dir / json_file.name
            if not dest_json.exists():
                shutil.copy(json_file, dest_json)

            # Copy network_setup.json into job dir
            dest_network = job_dir / "network_setup.json"
            if not dest_network.exists():
                shutil.copy(network_json, dest_network)

            # Create quickrun_output/ with a placeholder so HTCondor
            # always has a non-empty directory to transfer
            qo = job_dir / "quickrun_output"
            qo.mkdir(exist_ok=True)
            (qo / ".placeholder").touch()

            # Create empty result.json so transfer never fails
            (job_dir / "result.json").touch()

            # Copy the executable script into the job dir so it can be
            # listed in transfer_input_files (paths relative to initialdir).
            # This ensures the script reaches the execute node regardless
            # of how HTCondor resolves executable paths with initialdir.
            script_src = Path(args.script)
            dest_script = job_dir / script_src.name
            if not dest_script.exists():
                shutil.copy(script_src, dest_script)

            transform_entries.append(f"{cluster_id},{transform_name}")
            n_total += 1

    list_path = outdir / "transform_list.txt"
    with open(list_path, "w") as f:
        f.write("\n".join(transform_entries) + "\n")

    print(f"\nTotal transformation jobs: {n_total}")
    print(f"Wrote {list_path}")

    # Verify a sample job directory
    if transform_entries:
        parts = transform_entries[0].split(",")
        sample_dir = outdir / parts[0] / parts[1]
        print(f"\nSample job directory ({sample_dir}):")
        for p in sorted(sample_dir.rglob("*")):
            rel = p.relative_to(sample_dir)
            print(f"  {rel}")


if __name__ == "__main__":
    main()
