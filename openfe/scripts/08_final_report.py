#!/usr/bin/env python3
"""
openfe/scripts/08_final_report.py

Generate a comprehensive per-leg report for the full RBFE production
campaign. Covers all 1066 transformation legs with:

  - Status (COMPLETED / NaN / IndexError / OpenMMException / other / no_result)
  - Failure reason and a short message if failed
  - ddG estimate and uncertainty (for completed legs)
  - Wall-clock time (from HTCondor .log TimeExecute, or script timing)
  - GPU card type (from .log machine classad)
  - GPU capability and memory
  - Server / slot name (from .log SlotName)
  - Number of job starts / retries (from .log)
  - Result location (top-level result.json vs quickrun_output/result.json)

Also produces edge-level summary (complex+solvent pairing -> ddG) and
aggregate statistics by GPU type, server, and failure category.

Usage:
    python openfe/scripts/08_final_report.py \
        --production-dir openfe/production \
        --outdir openfe
"""

import argparse
import json
import re
from pathlib import Path

import pandas as pd


def load_result(result_path):
    """Return (estimate, uncertainty, status, fail_type, fail_msg)."""
    if not result_path.exists() or result_path.stat().st_size <= 100:
        return None, None, "NO_RESULT", None, None
    try:
        d = json.loads(result_path.read_text())
    except Exception as e:
        return None, None, "INVALID_JSON", "JSONDecodeError", str(e)[:120]

    est = d.get("estimate")
    unc = d.get("uncertainty")
    est_mag = est.get("magnitude") if est else None
    unc_mag = unc.get("magnitude") if unc else None

    if est_mag is not None:
        return est_mag, unc_mag, "COMPLETED", None, None

    # Null estimate - find the failure type from unit_results exceptions
    fail_type = "unknown"
    fail_msg = None
    for k, v in d.get("unit_results", {}).items():
        exc = v.get("exception")
        if exc and isinstance(exc, list) and len(exc) > 0:
            fail_type = exc[0]
            if len(exc) > 1:
                msg = exc[1]
                fail_msg = (msg[0] if isinstance(msg, list) and msg
                            else str(msg))[:120]
            break
    return None, None, "FAILED", fail_type, fail_msg


def parse_log(log_path):
    """Extract GPU/server/timing info from HTCondor .log file."""
    info = {
        "gpu_name": None, "gpu_capability": None, "gpu_memory_mb": None,
        "server": None, "slot": None, "time_execute_s": None,
        "num_job_starts": None, "exit_code": None,
    }
    if not log_path.exists():
        return info
    text = log_path.read_text(errors="replace")

    m = re.search(r'DeviceName = "([^"]+)"', text)
    info["gpu_name"] = m.group(1) if m else None
    m = re.search(r'Capability = ([\d.]+)', text)
    info["gpu_capability"] = float(m.group(1)) if m else None
    m = re.search(r'GlobalMemoryMb = (\d+)', text)
    info["gpu_memory_mb"] = int(m.group(1)) if m else None

    # SlotName: slot2_3@gitter0000.chtc.wisc.edu
    m = re.search(r'SlotName:\s*(\S+)', text)
    if m:
        info["slot"] = m.group(1)
        if "@" in m.group(1):
            info["server"] = m.group(1).split("@", 1)[1]

    # Last TimeExecute (final run if multiple)
    times = re.findall(r'TimeExecute \(s\)\s*:\s*(\d+)', text)
    info["time_execute_s"] = int(times[-1]) if times else None

    # Count job execution events (number of starts)
    info["num_job_starts"] = len(re.findall(r'Job executing on host', text))

    # Last exit code
    codes = re.findall(r'Normal termination \(return value (\d+)\)', text)
    info["exit_code"] = int(codes[-1]) if codes else None

    return info


def parse_out_timing(job_dir):
    """Get wall-clock and host from script .out if present."""
    wall_s = None
    host = None
    for out_file in job_dir.glob("**/*.out"):
        text = out_file.read_text(errors="replace")
        for line in text.splitlines():
            if "Wall-clock:" in line:
                try:
                    wall_s = int(line.split("Wall-clock:")[1].split("s")[0].strip())
                except (ValueError, IndexError):
                    pass
            if line.startswith("Host:") and host is None:
                host = line.split("Host:", 1)[1].strip()
    return wall_s, host


def edge_name(transform_name):
    """rbfe_<A>_complex_<B>_complex -> (rbfe_<A>_<B>, leg)."""
    if transform_name.endswith("_complex"):
        leg = "complex"
    elif transform_name.endswith("_solvent"):
        leg = "solvent"
    else:
        return transform_name, None
    edge = re.sub(r'_(complex|solvent)$', '', transform_name)
    edge = re.sub(r'_(complex|solvent)_', '_', edge)
    return edge, leg


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--production-dir", default="openfe/production")
    ap.add_argument("--outdir", default="openfe")
    args = ap.parse_args()

    prod_dir = Path(args.production_dir)
    outdir = Path(args.outdir)
    transform_list = prod_dir / "transform_list.txt"

    rows = []
    with open(transform_list) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            cluster_id, transform_name = [x.strip() for x in line.split(",", 1)]
            job_dir = prod_dir / cluster_id / transform_name

            # Try top-level result.json first, then quickrun_output/
            top = job_dir / "result.json"
            inner = job_dir / "quickrun_output" / "result.json"
            est, unc, status, ftype, fmsg = load_result(top)
            result_loc = "top"
            if status != "COMPLETED":
                est2, unc2, status2, ftype2, fmsg2 = load_result(inner)
                if status2 == "COMPLETED":
                    est, unc, status, ftype, fmsg = est2, unc2, status2, ftype2, fmsg2
                    result_loc = "inner"
                elif status == "NO_RESULT" and status2 != "NO_RESULT":
                    # inner has failure info, top doesn't
                    est, unc, status, ftype, fmsg = est2, unc2, status2, ftype2, fmsg2
                    result_loc = "inner"

            log_files = list(job_dir.glob("**/*.log"))
            log_info = parse_log(log_files[0]) if log_files else {}
            wall_s, host = parse_out_timing(job_dir)

            edge, leg = edge_name(transform_name)

            rows.append({
                "cluster_id": cluster_id,
                "transform_name": transform_name,
                "edge": edge,
                "leg": leg,
                "status": status,
                "fail_type": ftype,
                "fail_msg": fmsg,
                "ddg_estimate": est,
                "ddg_uncertainty": unc,
                "result_location": result_loc if status == "COMPLETED" else None,
                "wall_clock_s": wall_s or log_info.get("time_execute_s"),
                "gpu_name": log_info.get("gpu_name"),
                "gpu_capability": log_info.get("gpu_capability"),
                "gpu_memory_mb": log_info.get("gpu_memory_mb"),
                "server": log_info.get("server"),
                "slot": log_info.get("slot"),
                "num_job_starts": log_info.get("num_job_starts"),
                "exit_code": log_info.get("exit_code"),
            })

    df = pd.DataFrame(rows)
    df["wall_clock_h"] = df["wall_clock_s"] / 3600

    # Write the full per-leg report
    report_path = outdir / "final_leg_report.csv"
    df.to_csv(report_path, index=False)

    # ============ Console summary ============
    print("=" * 70)
    print("FINAL CAMPAIGN REPORT")
    print("=" * 70)
    print(f"Total legs: {len(df)}")
    print()
    print("STATUS:")
    print(df["status"].value_counts().to_string())

    print("\nFAILURE TYPES:")
    failed = df[df["status"] != "COMPLETED"]
    if len(failed):
        print(failed["fail_type"].value_counts().to_string())

    print("\nBY LEG TYPE:")
    for leg in ["complex", "solvent"]:
        sub = df[df["leg"] == leg]
        done = (sub["status"] == "COMPLETED").sum()
        print(f"  {leg}: {done}/{len(sub)} completed")

    # Timing
    completed = df[df["status"] == "COMPLETED"]
    if len(completed) and completed["wall_clock_h"].notna().any():
        wc = completed["wall_clock_h"].dropna()
        print(f"\nTIMING (completed, n={len(wc)}):")
        print(f"  mean={wc.mean():.1f}h  median={wc.median():.1f}h  "
              f"min={wc.min():.1f}h  max={wc.max():.1f}h")

    # By GPU type
    print("\nBY GPU TYPE:")
    gpu_stats = df.groupby("gpu_name").agg(
        n_legs=("transform_name", "count"),
        n_completed=("status", lambda x: (x == "COMPLETED").sum()),
        mean_h=("wall_clock_h", "mean"),
    ).reset_index().sort_values("n_legs", ascending=False)
    print(gpu_stats.to_string(index=False))

    # By server (top 15)
    print("\nBY SERVER (top 15 by leg count):")
    srv_stats = df.groupby("server").agg(
        n_legs=("transform_name", "count"),
        n_completed=("status", lambda x: (x == "COMPLETED").sum()),
        n_failed=("status", lambda x: (x != "COMPLETED").sum()),
    ).reset_index().sort_values("n_legs", ascending=False).head(15)
    print(srv_stats.to_string(index=False))

    # ============ Edge-level summary ============
    edge_rows = []
    for edge, grp in df.groupby(["cluster_id", "edge"]):
        complex_leg = grp[grp["leg"] == "complex"]
        solvent_leg = grp[grp["leg"] == "solvent"]
        c_done = len(complex_leg) and complex_leg.iloc[0]["status"] == "COMPLETED"
        s_done = len(solvent_leg) and solvent_leg.iloc[0]["status"] == "COMPLETED"
        ddg = None
        if c_done and s_done:
            ddg = (complex_leg.iloc[0]["ddg_estimate"]
                   - solvent_leg.iloc[0]["ddg_estimate"])
        edge_rows.append({
            "cluster_id": edge[0],
            "edge": edge[1],
            "complex_done": bool(c_done),
            "solvent_done": bool(s_done),
            "both_done": bool(c_done and s_done),
            "ddg": ddg,
            "complex_fail": (complex_leg.iloc[0]["fail_type"]
                             if len(complex_leg) and not c_done else None),
            "solvent_fail": (solvent_leg.iloc[0]["fail_type"]
                             if len(solvent_leg) and not s_done else None),
        })
    edge_df = pd.DataFrame(edge_rows)
    edge_path = outdir / "final_edge_report.csv"
    edge_df.to_csv(edge_path, index=False)

    print("\n" + "=" * 70)
    print("EDGE-LEVEL SUMMARY")
    print("=" * 70)
    print(f"Total edges: {len(edge_df)}")
    print(f"Both legs done (usable ddG): {edge_df['both_done'].sum()}")
    print(f"Incomplete edges: {(~edge_df['both_done']).sum()}")
    if edge_df["both_done"].any():
        ddgs = edge_df[edge_df["both_done"]]["ddg"].dropna()
        print(f"\nddG distribution (kcal/mol):")
        print(f"  mean={ddgs.mean():.2f}  median={ddgs.median():.2f}  "
              f"min={ddgs.min():.2f}  max={ddgs.max():.2f}")
        print(f"  |ddG|>5: {(ddgs.abs()>5).sum()}  "
              f"|ddG|>10: {(ddgs.abs()>10).sum()}")

    print(f"\nWrote per-leg report:  {report_path}")
    print(f"Wrote per-edge report: {edge_path}")


if __name__ == "__main__":
    main()
