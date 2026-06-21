#!/usr/bin/env python3
"""
openfe/scripts/06_monitor_production.py

Monitor the production openfe quickrun campaign and optionally clean
up large checkpoint files from completed jobs.

For each job directory, checks:
  - result.json size > 100 bytes  -> completed successfully
  - quickrun_output/ with checkpoints -> in progress or failed
  - Empty result.json + empty quickrun_output/ -> not started

Also reports storage usage and can clean up large files from completed
jobs that weren't cleaned by the shell script (e.g., if eviction
happened during cleanup).

Usage:
    python openfe/scripts/06_monitor_production.py \
        --production-dir openfe/production

    # With cleanup of completed job checkpoint files:
    python openfe/scripts/06_monitor_production.py \
        --production-dir openfe/production --cleanup
"""

import argparse
import json
import os
from pathlib import Path

import pandas as pd


def get_dir_size_mb(path):
    """Get total size of a directory in MB."""
    total = 0
    for dirpath, dirnames, filenames in os.walk(path):
        for f in filenames:
            fp = os.path.join(dirpath, f)
            try:
                total += os.path.getsize(fp)
            except OSError:
                pass
    return total / (1024 * 1024)


def check_job(job_dir):
    """Check completion status of a single job directory."""
    result_file = job_dir / "result.json"
    qo_dir = job_dir / "quickrun_output"

    status = "NOT_STARTED"
    wall_clock_s = None
    host = None
    size_mb = get_dir_size_mb(job_dir) if job_dir.exists() else 0

    if result_file.exists() and result_file.stat().st_size > 100:
        status = "COMPLETED"
        # Find .out file - may be condor.out or quickrun.<cluster>.<proc>.out
        out_files = list(job_dir.glob("**/*.out"))
        out_file = out_files[0] if out_files else None
        if out_file and out_file.exists():
            text = out_file.read_text(errors="replace")
            for line in text.splitlines():
                if "Wall-clock:" in line:
                    try:
                        wall_clock_s = int(line.split("Wall-clock:")[1].split("s")[0].strip())
                    except (ValueError, IndexError):
                        pass
                if "Host:" in line and host is None:
                    host = line.split("Host:")[1].strip()

        # Check if large files still need cleanup
        has_large_files = False
        if qo_dir.exists():
            for p in qo_dir.rglob("*.nc"):
                has_large_files = True
                break
            if not has_large_files:
                for p in qo_dir.rglob("*.chk"):
                    has_large_files = True
                    break
        if has_large_files:
            status = "COMPLETED_NEEDS_CLEANUP"

    elif qo_dir.exists():
        cache_dir = qo_dir / "quickrun_cache"
        has_sim = any(qo_dir.glob("shared_*SimulationUnit*"))
        has_setup = any(qo_dir.glob("shared_*SetupUnit*"))
        if cache_dir.exists() or has_sim or has_setup:
            status = "IN_PROGRESS_OR_FAILED"
        elif any(qo_dir.iterdir()):
            status = "STARTED"

    return {
        "status": status,
        "wall_clock_s": wall_clock_s,
        "wall_clock_h": wall_clock_s / 3600 if wall_clock_s else None,
        "host": host,
        "size_mb": round(size_mb, 1),
    }


def cleanup_completed(job_dir):
    """Delete large checkpoint/trajectory files from a completed job."""
    qo_dir = job_dir / "quickrun_output"
    if not qo_dir.exists():
        return 0
    freed = 0
    for pattern in ["shared_*SimulationUnit*", "shared_*SetupUnit*"]:
        for d in qo_dir.glob(pattern):
            if d.is_dir():
                size = get_dir_size_mb(d)
                import shutil
                shutil.rmtree(d)
                freed += size
    for ext in ["*.nc", "*.chk"]:
        for f in qo_dir.rglob(ext):
            freed += f.stat().st_size / (1024 * 1024)
            f.unlink()
    return freed


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--production-dir", default="openfe/production")
    ap.add_argument("--outdir", default="openfe")
    ap.add_argument("--cleanup", action="store_true",
                    help="Delete large checkpoint files from completed jobs")
    ap.add_argument("--cleanup-failed", action="store_true",
                    help="Delete large checkpoint files from permanently "
                         "failed jobs (those with IN_PROGRESS_OR_FAILED "
                         "status and no running HTCondor job). Use with "
                         "caution: only run after confirming these jobs "
                         "are truly done retrying (condor_q shows no "
                         "matching jobs).")
    args = ap.parse_args()

    prod_dir = Path(args.production_dir)
    outdir = Path(args.outdir)

    transform_list = prod_dir / "transform_list.txt"
    if not transform_list.exists():
        print(f"ERROR: {transform_list} not found")
        return

    rows = []
    with open(transform_list) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = [x.strip() for x in line.split(",", 1)]
            if len(parts) != 2:
                continue
            cluster_id, transform_name = parts[0], parts[1]
            job_dir = prod_dir / cluster_id / transform_name

            result = check_job(job_dir)
            result["cluster_id"] = cluster_id
            result["transform_name"] = transform_name
            result["is_complex"] = "_complex_" in transform_name
            rows.append(result)

    df = pd.DataFrame(rows)

    # Summary
    print("=" * 60)
    print("PRODUCTION CAMPAIGN STATUS")
    print("=" * 60)
    print(f"Total transformations: {len(df)}")
    print()
    print(df["status"].value_counts().to_string())

    n_completed = df["status"].isin(["COMPLETED", "COMPLETED_NEEDS_CLEANUP"]).sum()
    n_total = len(df)
    pct = 100 * n_completed / n_total if n_total > 0 else 0
    print(f"\nProgress: {n_completed} / {n_total} ({pct:.1f}%)")

    # By leg type
    print("\nBy leg type:")
    for is_complex in [True, False]:
        leg_name = "complex" if is_complex else "solvent"
        sub = df[df["is_complex"] == is_complex]
        n_done = sub["status"].isin(["COMPLETED", "COMPLETED_NEEDS_CLEANUP"]).sum()
        print(f"  {leg_name}: {n_done} / {len(sub)}")

    # Storage
    total_mb = df["size_mb"].sum()
    print(f"\nTotal storage used: {total_mb:.0f} MB ({total_mb/1024:.1f} GB)")

    # Timing for completed jobs
    completed = df[df["status"].isin(["COMPLETED", "COMPLETED_NEEDS_CLEANUP"])]
    if len(completed) and completed["wall_clock_h"].notna().any():
        wc = completed["wall_clock_h"].dropna()
        print(f"\nCompleted job timing (n={len(wc)}):")
        print(f"  Mean: {wc.mean():.1f}h")
        print(f"  Min:  {wc.min():.1f}h")
        print(f"  Max:  {wc.max():.1f}h")

        remaining = n_total - n_completed
        if remaining > 0:
            mean_h = wc.mean()
            print(f"\nRemaining: {remaining} jobs")
            for n_gpus in [10, 20, 30, 50]:
                est_days = (remaining * mean_h) / (n_gpus * 24)
                print(f"  At {n_gpus} GPUs: {est_days:.1f} days")

    # Per-cluster completion
    cluster_status = df.groupby("cluster_id").agg(
        n_total=("status", "count"),
        n_done=("status", lambda x: x.isin(["COMPLETED", "COMPLETED_NEEDS_CLEANUP"]).sum()),
    ).reset_index()
    cluster_status["pct"] = 100 * cluster_status["n_done"] / cluster_status["n_total"]

    fully_done = (cluster_status["n_done"] == cluster_status["n_total"]).sum()
    print(f"\nClusters fully completed: {fully_done} / {len(cluster_status)}")

    # Cleanup completed jobs
    if args.cleanup:
        needs_cleanup = df[df["status"] == "COMPLETED_NEEDS_CLEANUP"]
        if len(needs_cleanup):
            print(f"\nCleaning up {len(needs_cleanup)} completed jobs "
                  f"with remaining large files...")
            total_freed = 0
            for _, row in needs_cleanup.iterrows():
                job_dir = prod_dir / row["cluster_id"] / row["transform_name"]
                freed = cleanup_completed(job_dir)
                total_freed += freed
            print(f"Freed {total_freed:.0f} MB")
        else:
            print("\nNo completed jobs need cleanup.")

    # Cleanup permanently failed jobs
    if args.cleanup_failed:
        failed = df[df["status"] == "IN_PROGRESS_OR_FAILED"]
        if len(failed):
            print(f"\nWARNING: cleaning up {len(failed)} failed jobs' "
                  f"checkpoint files. These cannot be resumed after "
                  f"cleanup.")
            total_freed = 0
            for _, row in failed.iterrows():
                job_dir = prod_dir / row["cluster_id"] / row["transform_name"]
                freed = cleanup_completed(job_dir)
                total_freed += freed
            print(f"Freed {total_freed:.0f} MB")
        else:
            print("\nNo failed jobs to clean up.")

    # Generate retry list for incomplete jobs
    incomplete = df[~df["status"].isin(["COMPLETED", "COMPLETED_NEEDS_CLEANUP"])]
    if len(incomplete) > 0:
        retry_path = prod_dir / "transform_list_retry.txt"
        with open(retry_path, "w") as f:
            for _, row in incomplete.iterrows():
                f.write(f"{row['cluster_id']},{row['transform_name']}\n")
        print(f"\nWrote retry list: {retry_path} ({len(incomplete)} jobs)")
        print("To resubmit failed/incomplete jobs:")
        print(f"  cd openfe")
        print(f"  # First remove any held jobs from the current submission:")
        print(f"  condor_rm <cluster_id>")
        print(f"  # Then resubmit just the incomplete jobs:")
        print(f'  condor_submit submit_quickrun_production.sub \\')
        print(f'    -append "queue cluster_id,transform_name from '
              f'production/transform_list_retry.txt"')
        print(f"\n  # Or with server exclusions if specific servers are broken:")
        print(f'  condor_submit submit_quickrun_production.sub \\')
        print(f'    -append \'requirements = (Machine != '
              f'"bad-server.chtc.wisc.edu")\' \\')
        print(f'    -append "queue cluster_id,transform_name from '
              f'production/transform_list_retry.txt"')

    # Save status
    out_path = outdir / "production_status.csv"
    df.to_csv(out_path, index=False)
    print(f"\nWrote {out_path}")


if __name__ == "__main__":
    main()
