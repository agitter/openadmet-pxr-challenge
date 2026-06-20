#!/usr/bin/env python3
"""
openfe/scripts/analyze_quickrun_timing.py

Analyze openfe quickrun timing test results across two HTCondor
submissions (cluster IDs 7957251 and 7966096), using only the
log/out/err files which are reliably unique per job.

The timing_test/ output directory had collisions (second submission
overwrote first), so timing data from timing_test/ is used only as
secondary confirmation, matched by hostname to avoid ambiguity.

Per job, extracts from HTCondor files:
  .log  -> server hostname, GPU device name, capability, memory,
           job exit status, wall-clock time (TimeExecute)
  .out  -> CUDA available (True/False), iterations completed,
           estimated/actual completion time from openmmtools progress
  .err  -> failure reason (CUDA not found, NaN, etc.)

Outputs a summary table and per-GPU-type timing statistics.

Usage:
    python openfe/scripts/analyze_quickrun_timing.py \
        --log-dir openfe/logs \
        --timing-dir openfe/timing_test \
        --clusters 7957251 7966096
"""

import argparse
import re
from pathlib import Path

import pandas as pd


# ---------------------------------------------------------------
# Parsers for each file type
# ---------------------------------------------------------------

def parse_log(path):
    """Parse HTCondor .log file for GPU info, exit status, wall-clock."""
    text = path.read_text(errors="replace")
    result = {}

    # GPU device info from the machine classad in job executing event
    m = re.search(r'DeviceName = "([^"]+)"', text)
    result["gpu_name"] = m.group(1) if m else None

    m = re.search(r'Capability = ([\d.]+)', text)
    result["gpu_capability"] = float(m.group(1)) if m else None

    m = re.search(r'GlobalMemoryMb = (\d+)', text)
    result["gpu_memory_mb"] = int(m.group(1)) if m else None

    # Execute host from SlotName line (most reliable)
    # e.g. "SlotName: slot2_7@gpu2011.chtc.wisc.edu"
    m = re.search(r'SlotName: \S+@([\w.]+)', text)
    result["server"] = m.group(1) if m else None

    # Job exit status
    m = re.search(r'Normal termination \(return value (\d+)\)', text)
    result["exit_code"] = int(m.group(1)) if m else None

    # TimeExecute from resource usage summary
    m = re.search(r'TimeExecute \(s\)\s*:\s*(\d+)', text)
    result["time_execute_s"] = int(m.group(1)) if m else None

    return result


def parse_out(path):
    """Parse job .out file for CUDA availability and wall-clock time."""
    text = path.read_text(errors="replace")
    result = {}

    # CUDA check result from our fast-fail check
    if "CUDA platform confirmed available" in text:
        result["cuda_available"] = True
    elif "CUDA platform not available" in text:
        result["cuda_available"] = False
    else:
        result["cuda_available"] = None

    # Short hostname from "Host: <name>" line
    m = re.search(r'^Host: ([\w.]+)', text, re.MULTILINE)
    result["host_from_out"] = m.group(1) if m else None

    # Wall-clock time from our timing summary
    m = re.search(r'Wall-clock time: (\d+) seconds', text)
    result["wall_clock_s"] = int(m.group(1)) if m else None

    # Result dG if simulation completed
    m = re.search(r'dG = ([\d.+-]+) kilocalorie_per_mole', text)
    result["dG_kcal_mol"] = float(m.group(1)) if m else None

    return result


def parse_err(path):
    """Parse job .err file for failure reason and iteration timing.
    
    Note: openfe writes 'CUDA-based GPU not found' as an informational
    probe message even on successful GPU runs - this is NOT a failure.
    Only treat it as a failure if the job also has NaN errors or our
    fast-fail CUDA check fired.
    """
    text = path.read_text(errors="replace")
    result = {}

    # Hostname from openfe system probe (more reliable than .log alias)
    m = re.search(r"hostname: '([\w.]+)'", text)
    result["host_from_err"] = m.group(1) if m else None

    # Real failure indicators (not just the probe message)
    has_nan = "SimulationNaNError" in text or "resulted in a NaN" in text
    has_fast_fail = "CUDA platform not available" in text  # our check
    has_protocol_error = "The protocol unit" in text and "failed" in text

    if has_fast_fail:
        result["failure_reason"] = "CUDA not available (fast-fail)"
    elif has_nan:
        result["failure_reason"] = "NaN error (ran on CPU without CUDA)"
    elif has_protocol_error:
        m = re.search(r'Error: (.+)', text)
        result["failure_reason"] = m.group(1)[:80] if m else "Protocol error"
    else:
        result["failure_reason"] = None

    # Iteration timing - openmmtools logs go to stderr
    iter_times = re.findall(r'Iteration took ([\d.]+)s', text)
    if iter_times:
        times = [float(t) for t in iter_times]
        result["mean_s_per_iter"] = sum(times) / len(times)
        result["n_iter_timed"] = len(times)
    else:
        result["mean_s_per_iter"] = None
        result["n_iter_timed"] = None

    # Max iteration reached (also in stderr)
    iterations = re.findall(r'Iteration (\d+)/2000', text)
    result["max_iteration_err"] = int(max(iterations, key=int)) if iterations else None

    return result


# ---------------------------------------------------------------
# Main
# ---------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--log-dir", default="openfe/logs")
    ap.add_argument("--timing-dir", default="openfe/timing_test")
    ap.add_argument("--clusters", nargs="+", default=["7957251", "7966096"])
    ap.add_argument("--outdir", default="openfe")
    args = ap.parse_args()

    log_dir = Path(args.log_dir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    rows = []

    for cluster in args.clusters:
        # Find all process IDs for this cluster
        log_files = sorted(log_dir.glob(f"timing.{cluster}.*.log"))
        print(f"\nCluster {cluster}: {len(log_files)} jobs")

        for log_path in log_files:
            proc = re.search(r'\.(\d+)\.log$', log_path.name).group(1)
            out_path = log_dir / f"timing.{cluster}.{proc}.out"
            err_path = log_dir / f"timing.{cluster}.{proc}.err"

            row = {"cluster_id": cluster, "proc_id": proc,
                   "job_id": f"{cluster}.{proc}"}

            log_data = parse_log(log_path) if log_path.exists() else {}
            out_data = parse_out(out_path) if out_path.exists() else {}
            err_data = parse_err(err_path) if err_path.exists() else {}

            row.update(log_data)
            row.update(out_data)
            row.update(err_data)

            # wall_clock_s: from our timing summary in .out (most accurate)
            # fall back to TimeExecute from .log
            if not row.get("wall_clock_s") and row.get("time_execute_s"):
                row["wall_clock_s"] = row["time_execute_s"]

            # Use max_iteration from .err (openmmtools logs to stderr)
            if row.get("max_iteration_err"):
                row["max_iteration"] = row["max_iteration_err"]

            # Determine overall status using .out as ground truth.
            # The .err file may contain NaN errors from earlier simulation
            # attempts that OpenFE retried automatically - presence of a
            # dG result in .out means the job ultimately succeeded.
            completed = row.get("dG_kcal_mol") is not None
            if completed:
                row["status"] = "SUCCESS"
                row["failure_reason"] = None  # clear any spurious failure
            elif row.get("failure_reason"):
                row["status"] = "FAILED"
            elif row.get("exit_code") == 0:
                # Exited 0 but no dG result - unusual, flag as unknown
                row["status"] = "UNKNOWN"
            elif row.get("exit_code") is not None:
                row["status"] = "FAILED"
            else:
                row["status"] = "UNKNOWN"

            # Projected total time: for completed jobs use actual wall_clock_s.
            # For jobs that ran but didn't finish, use mean_s_per_iter * 2000.
            if row.get("wall_clock_s") and completed:
                row["projected_total_s"] = row["wall_clock_s"]
                row["projected_total_h"] = row["wall_clock_s"] / 3600
            elif row.get("mean_s_per_iter"):
                row["projected_total_s"] = row["mean_s_per_iter"] * 2000
                row["projected_total_h"] = row["projected_total_s"] / 3600
            else:
                row["projected_total_s"] = None
                row["projected_total_h"] = None

            rows.append(row)

            # Print per-job summary
            status = row["status"]
            gpu = row.get("gpu_name") or "unknown GPU"
            server = row.get("server") or "unknown server"
            reason = row.get("failure_reason") or ""
            wc = f"{row['wall_clock_s']}s" if row.get("wall_clock_s") else "n/a"
            proj = f"{row['projected_total_h']:.1f}h" if row.get("projected_total_h") else "n/a"
            print(f"  {row['job_id']:15s} {status:8s} {gpu:35s} {server:30s} "
                  f"wall={wc:8s} proj={proj:6s} {reason}")

    df = pd.DataFrame(rows)

    # ---------------------------------------------------------------
    # Summary tables
    # ---------------------------------------------------------------
    print("\n" + "="*70)
    print("STATUS SUMMARY")
    print("="*70)
    print(df["status"].value_counts().to_string())

    print("\n" + "="*70)
    print("CUDA AVAILABILITY BY SERVER")
    print("="*70)
    df["cuda_worked"] = (df["status"] == "SUCCESS") & df["failure_reason"].isna()
    server_cuda = df.groupby(["server", "gpu_name", "gpu_capability"]).agg(
        cuda_worked=("cuda_worked", "sum"),
        total=("job_id", "count")
    ).reset_index()
    print(server_cuda.sort_values("gpu_capability").to_string(index=False))

    print("\n" + "="*70)
    print("TIMING BY GPU TYPE (successful jobs only)")
    print("="*70)
    success = df[df["status"] == "SUCCESS"].copy()
    if len(success):
        by_gpu = success.groupby(["gpu_name","gpu_capability"]).agg(
            n_jobs=("job_id","count"),
            mean_projected_h=("projected_total_h","mean"),
            min_projected_h=("projected_total_h","min"),
            max_projected_h=("projected_total_h","max"),
            mean_s_per_iter=("mean_s_per_iter","mean"),
        ).reset_index().sort_values("gpu_capability")
        print(by_gpu.to_string(index=False))

        print(f"\nOverall mean projected leg time: "
              f"{success['projected_total_h'].mean():.1f}h")
        print(f"Overall min: {success['projected_total_h'].min():.1f}h")
        print(f"Overall max: {success['projected_total_h'].max():.1f}h")

        print("\n" + "="*70)
        print("PRODUCTION CAMPAIGN ESTIMATE (1066 total legs)")
        print("="*70)
        mean_h = success["projected_total_h"].mean()
        for n_gpus in [10, 20, 30, 50]:
            days = (1066 * mean_h) / (n_gpus * 24)
            print(f"  {n_gpus:3d} GPUs: {days:.1f} days")
    else:
        print("No successful jobs to report timing for.")

    print("\n" + "="*70)
    print("FAILED JOBS")
    print("="*70)
    failed = df[df["status"] == "FAILED"]
    if len(failed):
        print(failed[["job_id","server","gpu_name","gpu_capability",
                       "failure_reason"]].to_string(index=False))
    else:
        print("No failed jobs.")

    # Save full table
    out_path = outdir / "timing_analysis.csv"
    df.to_csv(out_path, index=False)
    print(f"\nWrote full table -> {out_path}")


if __name__ == "__main__":
    main()
