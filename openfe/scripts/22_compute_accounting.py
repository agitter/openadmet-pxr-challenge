#!/usr/bin/env python3
"""
openfe/scripts/22_compute_accounting.py

Comprehensive GPU compute accounting for the RBFE campaign. Parses every
HTCondor event log (.log) under the production tree - all attempts of all
legs, completed and failed, Kartograf and salvage - and reports the true
total GPU burn.

HTCondor event-log format (NOT the resource-usage summary):
  001 ... Job executing on host: ...   <- attempt start (+ timestamp)
         SlotName: <prefix>@<host>     <- slot identity
         GPUs_GPU_xxx = [ ...; Capability = X; DeviceName = "..."; ... ]
  004 ... Job was evicted.             <- attempt end (+ timestamp)
         Run Remote Usage  Usr ... Sys ...
  005 ... Job terminated.              <- final attempt end
Each attempt's wall time = (end_ts - start_ts). Total burn sums ALL
attempts (the many-attempt churn from CUDA fast-fail evictions is real
cost). One leg = one transform dir; its attempts may span several logs
(multiple submission clusters).

Slot classification (SlotName prefix authoritative, then host):
  glidein_*   -> OSG/OSPool
  backfill*   -> backfill
  slotN_M @ {gitter0000,gitter2003,gpulab2001} -> prioritized (BMI_Gitter)
  slotN_M @ other *.chtc.wisc.edu / researcher nodes -> shared (CHTC GPU Lab)

GPU grouped TWO independent ways: by DeviceName, and by Capability.

Usage:
    python openfe/scripts/22_compute_accounting.py \
        --production-dir openfe/production \
        --outdir openfe
"""

import argparse
import re
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PRIORITIZED_HOSTS = {"gitter0000.chtc.wisc.edu", "gitter2003.chtc.wisc.edu",
                     "gpulab2001.chtc.wisc.edu"}

# Event-log line patterns
RE_EXEC = re.compile(r'^001 \(([\d.]+)\) (\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}) '
                     r'Job executing on host')
RE_END = re.compile(r'^00[45] \(([\d.]+)\) (\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}) '
                    r'Job (was evicted|terminated)')
RE_SLOT = re.compile(r'SlotName:\s*(\S+)')
RE_GPU = re.compile(r'Capability = ([\d.]+);\s*CoresPerCU.*?'
                    r'DeviceName = "([^"]+)"')
# DeviceName may precede or follow Capability depending on attr order;
# use a more tolerant two-field search within a GPUs_ classad line.
RE_CAP = re.compile(r'Capability = ([\d.]+)')
RE_DEV = re.compile(r'DeviceName = "([^"]+)"')
RE_RETVAL = re.compile(r'Normal termination \(return value (\d+)\)')


def parse_ts(s):
    return datetime.strptime(s, "%Y-%m-%d %H:%M:%S")


def classify_slot(slotname):
    """Return one of: osg, backfill, prioritized, shared, unknown."""
    if slotname is None:
        return "unknown"
    prefix = slotname.split("@")[0]
    if prefix.startswith("glidein_") or "glidein_" in slotname:
        return "osg"
    if prefix.startswith("backfill"):
        return "backfill"
    # regular slotN_M -> test host
    host = slotname.split("@")[-1]
    if host in PRIORITIZED_HOSTS:
        return "prioritized"
    return "shared"


def parse_log_attempts(path):
    """Parse one .log into a list of attempt dicts, each with start/end
    timestamps, slot, host, GPU device + capability, and end reason."""
    lines = path.read_text(errors="replace").splitlines()
    attempts = []
    cur = None
    i = 0
    while i < len(lines):
        line = lines[i]
        m = RE_EXEC.match(line)
        if m:
            # close any open attempt without explicit end (shouldn't happen)
            cur = {"start": parse_ts(m.group(2)), "end": None,
                   "slot": None, "host": None, "device": None,
                   "capability": None, "end_reason": None, "retval": None}
            # scan following indented block for SlotName + GPU classad
            j = i + 1
            block = []
            while j < len(lines) and not re.match(r'^\d{3} \(', lines[j]):
                block.append(lines[j])
                j += 1
            btext = "\n".join(block)
            ms = RE_SLOT.search(btext)
            if ms:
                cur["slot"] = ms.group(1)
                cur["host"] = ms.group(1).split("@")[-1]
            mc = RE_CAP.search(btext)
            if mc:
                cur["capability"] = float(mc.group(1))
            md = RE_DEV.search(btext)
            if md:
                cur["device"] = md.group(1)
            attempts.append(cur)
            i = j
            continue
        me = RE_END.match(line)
        if me and attempts:
            # attach end to the most recent attempt lacking an end
            for a in reversed(attempts):
                if a["end"] is None:
                    a["end"] = parse_ts(me.group(2))
                    a["end_reason"] = me.group(3)
                    # look ahead for return value in the block
                    j = i + 1
                    blk = []
                    while j < len(lines) and not re.match(r'^\d{3} \(', lines[j]):
                        blk.append(lines[j]); j += 1
                    mr = RE_RETVAL.search("\n".join(blk))
                    if mr:
                        a["retval"] = int(mr.group(1))
                    break
        i += 1
    return attempts


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--production-dir", default="openfe/production")
    ap.add_argument("--outdir", default="openfe")
    args = ap.parse_args()

    prod = Path(args.production_dir)
    outdir = Path(args.outdir)

    log_files = sorted(prod.rglob("*.log"))
    print(f"Found {len(log_files)} .log files")

    rows = []
    for lp in log_files:
        # transform dir = parent of logs/ ; leg name = that dir's name
        # path: .../production/<cluster>/<transform>/logs/quickrun.X.Y.log
        transform_dir = lp.parent.parent
        leg = transform_dir.name
        is_salvage = "salvage" in str(lp).lower()
        leg_type = ("complex" if leg.endswith("_complex")
                    else "solvent" if leg.endswith("_solvent") else "unknown")
        for a in parse_log_attempts(lp):
            if a["end"] is not None:
                dur = (a["end"] - a["start"]).total_seconds()
            else:
                dur = None
            rows.append({
                "leg": leg, "leg_type": leg_type,
                "campaign": "salvage" if is_salvage else "kartograf",
                "log": lp.name, "host": a["host"], "slot": a["slot"],
                "slot_type": classify_slot(a["slot"]),
                "device": a["device"], "capability": a["capability"],
                "end_reason": a["end_reason"], "retval": a["retval"],
                "duration_s": dur,
            })

    df = pd.DataFrame(rows)
    df["duration_h"] = df["duration_s"] / 3600.0
    df.to_csv(outdir / "compute_accounting_attempts.csv", index=False)
    print(f"Parsed {len(df)} execution attempts across "
          f"{df['leg'].nunique()} legs")

    # ---- Total burn ----
    total_h = df["duration_h"].sum()
    print("\n" + "=" * 64)
    print("TOTAL GPU BURN (all attempts)")
    print("=" * 64)
    print(f"  Total GPU-hours (all attempts):  {total_h:,.1f}")
    print(f"  Total attempts:                  {len(df):,}")
    print(f"  Total legs:                      {df['leg'].nunique():,}")
    print(f"  Mean attempts per leg:           "
          f"{len(df)/df['leg'].nunique():.1f}")

    # final-attempt-only burn (last attempt per leg by start time proxy:
    # here, the longest or the terminated one). Approx: per leg, the attempt
    # whose end_reason=='terminated' is the final; else the longest.
    def final_attempt_h(g):
        term = g[g["end_reason"] == "terminated"]
        if len(term):
            return term["duration_h"].sum()
        return g["duration_h"].max() if len(g) else 0
    final_h = df.groupby("leg").apply(final_attempt_h).sum()
    print(f"\n  Final-attempt-only GPU-hours:    {final_h:,.1f}")
    print(f"  Retry/churn overhead:            {total_h - final_h:,.1f} "
          f"({100*(total_h-final_h)/total_h:.1f}% of total)")

    # ---- Burn by campaign ----
    print("\n" + "=" * 64)
    print("BURN BY CAMPAIGN")
    print("=" * 64)
    camp = df.groupby("campaign").agg(
        attempts=("duration_h", "size"),
        gpu_hours=("duration_h", "sum")).reset_index()
    print(camp.to_string(index=False))

    # ---- Burn by slot type ----
    print("\n" + "=" * 64)
    print("BURN BY SLOT TYPE")
    print("=" * 64)
    slot = df.groupby("slot_type").agg(
        attempts=("duration_h", "size"),
        gpu_hours=("duration_h", "sum"),
        mean_attempt_h=("duration_h", "mean")).reset_index()
    slot = slot.sort_values("gpu_hours", ascending=False)
    print(slot.to_string(index=False))

    # ---- Burn by specific GPU device ----
    print("\n" + "=" * 64)
    print("BURN BY GPU DEVICE TYPE")
    print("=" * 64)
    dev = df[df["device"].notna()].groupby("device").agg(
        attempts=("duration_h", "size"),
        gpu_hours=("duration_h", "sum"),
        mean_attempt_h=("duration_h", "mean")).reset_index()
    dev = dev.sort_values("gpu_hours", ascending=False)
    print(dev.to_string(index=False))

    # ---- Burn by capability ----
    print("\n" + "=" * 64)
    print("BURN BY GPU CAPABILITY")
    print("=" * 64)
    cap = df[df["capability"].notna()].groupby("capability").agg(
        attempts=("duration_h", "size"),
        gpu_hours=("duration_h", "sum"),
        mean_attempt_h=("duration_h", "mean")).reset_index()
    cap = cap.sort_values("capability")
    print(cap.to_string(index=False))

    # ---- Successful-leg timing ----
    # A "successful" attempt: terminated normally with retval 0, OR the
    # longest attempt of a leg that ultimately produced output. Here we use
    # attempts that ran a meaningful duration and terminated (not evicted).
    success = df[(df["end_reason"] == "terminated") &
                 (df["retval"] == 0)].copy()
    print("\n" + "=" * 64)
    print("SUCCESSFUL-LEG TIMING (terminated, retval 0)")
    print("=" * 64)
    print(f"  Successful terminations: {len(success)}")
    if len(success):
        print(f"  Mean: {success['duration_h'].mean():.2f}h  "
              f"median: {success['duration_h'].median():.2f}h  "
              f"min: {success['duration_h'].min():.2f}h  "
              f"max: {success['duration_h'].max():.2f}h")
        print("\n  By GPU device:")
        sd = success.groupby("device").agg(
            n=("duration_h", "size"), mean_h=("duration_h", "mean"),
            median_h=("duration_h", "median"),
            min_h=("duration_h", "min"),
            max_h=("duration_h", "max")).reset_index()
        print(sd.sort_values("mean_h").to_string(index=False))
        print("\n  By capability:")
        sc = success.groupby("capability").agg(
            n=("duration_h", "size"), mean_h=("duration_h", "mean"),
            median_h=("duration_h", "median")).reset_index()
        print(sc.sort_values("capability").to_string(index=False))

    # ---- Where jobs ran: by host ----
    print("\n" + "=" * 64)
    print("WHERE JOBS RAN (by host, top 25 by GPU-hours)")
    print("=" * 64)
    host = df[df["host"].notna()].groupby("host").agg(
        attempts=("duration_h", "size"),
        gpu_hours=("duration_h", "sum")).reset_index()
    host = host.sort_values("gpu_hours", ascending=False)
    print(host.head(25).to_string(index=False))

    # ---- Cost per usable result ----
    n_usable_edges = 316
    n_connected = 292
    print("\n" + "=" * 64)
    print("COST PER USABLE RESULT")
    print("=" * 64)
    print(f"  Total GPU-hours / usable edge (316):       "
          f"{total_h/n_usable_edges:.1f}")
    print(f"  Total GPU-hours / connected compound (292): "
          f"{total_h/n_connected:.1f}")

    # ---- Visualizations ----
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))

    # by slot type
    ax = axes[0, 0]
    ax.barh(slot["slot_type"], slot["gpu_hours"], color="#6baed6",
            edgecolor="k")
    ax.set_xlabel("GPU-hours"); ax.set_title("Burn by slot type")
    ax.invert_yaxis()

    # by device
    ax = axes[0, 1]
    ax.barh(dev["device"], dev["gpu_hours"], color="#74c476", edgecolor="k")
    ax.set_xlabel("GPU-hours"); ax.set_title("Burn by GPU device")
    ax.invert_yaxis()
    ax.tick_params(axis="y", labelsize=7)

    # by capability
    ax = axes[0, 2]
    ax.bar(cap["capability"].astype(str), cap["gpu_hours"],
           color="#fd8d3c", edgecolor="k")
    ax.set_xlabel("Capability"); ax.set_ylabel("GPU-hours")
    ax.set_title("Burn by capability")

    # attempts per leg distribution
    ax = axes[1, 0]
    apl = df.groupby("leg").size()
    ax.hist(apl, bins=range(1, apl.max() + 2), color="#9e9ac8",
            edgecolor="k", align="left")
    ax.set_xlabel("Attempts per leg"); ax.set_ylabel("Legs")
    ax.set_title(f"Multi-attempt churn (mean {apl.mean():.1f}/leg)")

    # successful timing by capability (box-ish: mean bars)
    ax = axes[1, 1]
    if len(success):
        sc2 = success.groupby("capability")["duration_h"].mean()
        ax.bar(sc2.index.astype(str), sc2.values, color="#41b6c4",
               edgecolor="k")
        ax.set_xlabel("Capability"); ax.set_ylabel("Mean successful leg (h)")
        ax.set_title("Successful-leg time by capability")

    # campaign split
    ax = axes[1, 2]
    ax.bar(camp["campaign"], camp["gpu_hours"], color=["#3182bd", "#e6550d"],
           edgecolor="k")
    ax.set_ylabel("GPU-hours"); ax.set_title("Burn: Kartograf vs salvage")

    plt.tight_layout()
    fig_path = outdir / "compute_accounting.png"
    plt.savefig(fig_path, dpi=130, bbox_inches="tight")
    print(f"\nWrote {fig_path}")
    print(f"Wrote {outdir/'compute_accounting_attempts.csv'}")


if __name__ == "__main__":
    main()
