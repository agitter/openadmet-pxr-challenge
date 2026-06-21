#!/usr/bin/env python3
"""
openfe/scripts/rescue_results.py

For jobs with IN_PROGRESS_OR_FAILED status, check whether
quickrun_output/result.json already contains a valid result
(from a completed first run where the shell script failed
to copy it due to the openfe -o file conflict error).

If found, copies quickrun_output/result.json to the top-level
result.json, effectively marking the job as completed without
needing to rerun anything.

Usage:
    python openfe/scripts/rescue_results.py \
        --production-dir openfe/production
"""

import argparse
import json
import shutil
from pathlib import Path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--production-dir", default="openfe/production")
    ap.add_argument("--dry-run", action="store_true",
                    help="Show what would be rescued without copying")
    args = ap.parse_args()

    prod_dir = Path(args.production_dir)
    transform_list = prod_dir / "transform_list.txt"

    n_rescued = 0
    n_checked = 0
    n_truly_failed = 0

    with open(transform_list) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            cluster_id, transform_name = [x.strip() for x in line.split(",", 1)]
            job_dir = prod_dir / cluster_id / transform_name

            top_result = job_dir / "result.json"
            inner_result = job_dir / "quickrun_output" / "result.json"

            # Only check jobs where top-level result.json is empty/missing
            if top_result.exists() and top_result.stat().st_size > 100:
                continue  # already completed

            n_checked += 1

            if inner_result.exists() and inner_result.stat().st_size > 100:
                # Validate it's real JSON with content
                try:
                    data = json.loads(inner_result.read_text())
                    if data:
                        if args.dry_run:
                            print(f"WOULD RESCUE: {cluster_id}/{transform_name}")
                        else:
                            shutil.copy(inner_result, top_result)
                            print(f"RESCUED: {cluster_id}/{transform_name}")
                        n_rescued += 1
                except (json.JSONDecodeError, Exception):
                    n_truly_failed += 1
            else:
                n_truly_failed += 1

    print(f"\nChecked {n_checked} incomplete jobs")
    print(f"Rescued: {n_rescued}")
    print(f"Truly failed (no inner result.json): {n_truly_failed}")
    if args.dry_run and n_rescued > 0:
        print("\nRe-run without --dry-run to perform the rescue")


if __name__ == "__main__":
    main()
