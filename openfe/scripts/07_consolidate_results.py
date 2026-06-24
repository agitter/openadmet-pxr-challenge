#!/usr/bin/env python3
"""
openfe/scripts/07_consolidate_results.py

Consolidate all RBFE leg results into a single consistent location and
verify completeness.

Background: results ended up in two places due to the campaign history:
  - Original 668 jobs (checkpoint-based run): quickrun_output/result.json
  - Rerun 398 jobs (no-checkpoint run): top-level result.json

This script:
  1. For each job, finds a valid result (non-null estimate) in either
     location, preferring the top-level result.json if both exist.
  2. Copies any valid result found only in quickrun_output/ up to the
     top-level result.json so all 1066 live in one place.
  3. Reports which jobs still have null/missing results (need rerun).
  4. Optionally removes quickrun_output/ directories once consolidated.

Usage:
    python openfe/scripts/07_consolidate_results.py \
        --production-dir openfe/production

    # After confirming all valid, clean up quickrun_output dirs:
    python openfe/scripts/07_consolidate_results.py \
        --production-dir openfe/production --cleanup
"""

import argparse
import json
import shutil
from pathlib import Path


def load_estimate(result_path):
    """Return the estimate magnitude (float) or None if null/invalid."""
    if not result_path.exists() or result_path.stat().st_size <= 100:
        return None
    try:
        d = json.loads(result_path.read_text())
        est = d.get("estimate")
        if est is not None and est.get("magnitude") is not None:
            return est["magnitude"]
    except Exception:
        pass
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--production-dir", default="openfe/production")
    ap.add_argument("--cleanup", action="store_true",
                    help="Remove quickrun_output/ directories after "
                         "consolidating valid results to top-level")
    args = ap.parse_args()

    prod_dir = Path(args.production_dir)
    transform_list = prod_dir / "transform_list.txt"
    if not transform_list.exists():
        print(f"ERROR: {transform_list} not found")
        return

    n_toplevel = 0       # valid result already at top level
    n_consolidated = 0   # valid result moved from quickrun_output to top
    n_null = 0           # no valid result anywhere
    null_jobs = []

    with open(transform_list) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            cluster_id, transform_name = [x.strip() for x in line.split(",", 1)]
            job_dir = prod_dir / cluster_id / transform_name

            top_result = job_dir / "result.json"
            inner_result = job_dir / "quickrun_output" / "result.json"

            top_est = load_estimate(top_result)
            inner_est = load_estimate(inner_result)

            if top_est is not None:
                # Already have a valid top-level result
                n_toplevel += 1
            elif inner_est is not None:
                # Valid result only in quickrun_output - copy it up
                shutil.copy(inner_result, top_result)
                n_consolidated += 1
            else:
                # No valid result anywhere
                n_null += 1
                null_jobs.append(f"{cluster_id},{transform_name}")

    total = n_toplevel + n_consolidated + n_null
    print("=" * 60)
    print("RESULT CONSOLIDATION")
    print("=" * 60)
    print(f"Total jobs:                  {total}")
    print(f"Already at top level:        {n_toplevel}")
    print(f"Consolidated from inner:     {n_consolidated}")
    print(f"Valid total:                 {n_toplevel + n_consolidated}")
    print(f"Null/missing (need rerun):   {n_null}")

    if null_jobs:
        retry_path = prod_dir / "transform_list_retry.txt"
        with open(retry_path, "w") as f:
            f.write("\n".join(null_jobs) + "\n")
        print(f"\nWrote {len(null_jobs)} null jobs to {retry_path}")
        if len(null_jobs) <= 20:
            for j in null_jobs:
                print(f"  {j}")
    else:
        print("\nAll 1066 jobs have valid results.")

    # Cleanup quickrun_output dirs
    if args.cleanup:
        if n_null > 0:
            print(f"\nNOT cleaning up: {n_null} jobs still need rerun. "
                  f"Their quickrun_output/ is preserved.")
        else:
            print("\nCleaning up quickrun_output/ directories...")
            freed_mb = 0
            n_removed = 0
            with open(transform_list) as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    cluster_id, transform_name = [
                        x.strip() for x in line.split(",", 1)]
                    qo = (prod_dir / cluster_id / transform_name
                          / "quickrun_output")
                    if qo.exists():
                        size = sum(p.stat().st_size
                                   for p in qo.rglob("*") if p.is_file())
                        freed_mb += size / (1024 * 1024)
                        shutil.rmtree(qo)
                        n_removed += 1
            print(f"Removed {n_removed} quickrun_output/ dirs, "
                  f"freed {freed_mb:.0f} MB")


if __name__ == "__main__":
    main()
