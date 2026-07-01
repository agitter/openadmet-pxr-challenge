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
# dsigpu* machines are also prioritized (a group the user belongs to)
PRIORITIZED_HOST_PREFIXES = ("dsigpu",)

# Threshold separating a genuine run attempt from the -o output-conflict
# requeue thrash (the dominant churn cause; see writeup). Attempts shorter
# than this that evicted are resume/output-conflict failures, not real work.
IMMEDIATE_FAIL_H = 5.0 / 60.0  # 5 minutes

# Event-log line patterns
RE_EXEC = re.compile(r'^001 \(([\d.]+)\) (\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}) '
                     r'Job executing on host')
RE_END = re.compile(r'^00[45] \(([\d.]+)\) (\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}) '
                    r'Job (was evicted|terminated)')
RE_SLOT = re.compile(r'SlotName:\s*(\S+)')
RE_CAP = re.compile(r'Capability = ([\d.]+)')
RE_DEV = re.compile(r'DeviceName = "([^"]+)"')
RE_RETVAL = re.compile(r'Normal termination \(return value (\d+)\)')


def parse_ts(s):
    return datetime.strptime(s, "%Y-%m-%d %H:%M:%S")


def classify_slot(slotname):
    """Return one of: osg, backfill, prioritized, shared, unknown.

    Rules (in order):
      backfill* prefix                      -> backfill
      glidein_* in name                     -> osg
      regular slotN_M on a non-CHTC host    -> osg  (OSPool sites:
          *-EP.* endpoints, montana.edu/nd.edu/amnh.org, short-names)
      regular slotN_M on dsigpu*/gitter0000/gitter2003/gpulab2001 (CHTC)
                                            -> prioritized
      regular slotN_M on other *.chtc.wisc.edu -> shared
    """
    if slotname is None:
        return "unknown"
    prefix = slotname.split("@")[0]
    host = slotname.split("@")[-1]
    if prefix.startswith("backfill"):
        return "backfill"
    if "glidein_" in slotname:
        return "osg"
    # Off-CHTC hosts reached with a regular slot are OSPool/OSG sites:
    #   OSPool endpoints (e.g. Colgate-CCARE-EP.*, GSU-Adonis-EP.*),
    #   off-site domains (montana.edu, nd.edu, amnh.org),
    #   and bare short-names (gpu03, rails04, compute-gpu-03, ...).
    if not host.endswith(".chtc.wisc.edu"):
        return "osg"
    # On-CHTC regular slots: prioritized hosts vs shared GPU Lab
    if host in PRIORITIZED_HOSTS or host.split(".")[0].startswith(
            PRIORITIZED_HOST_PREFIXES):
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
        # A leg "succeeded" if its transform dir holds a result.json
        leg_succeeded = (transform_dir / "result.json").exists()
        for a in parse_log_attempts(lp):
            if a["end"] is not None:
                dur = (a["end"] - a["start"]).total_seconds()
            else:
                dur = None
            dur_h = dur / 3600.0 if dur is not None else None
            # Tag attempt kind. immediate_fail = a short attempt that either
            # (a) was evicted (preemption) or (b) terminated with a nonzero
            # return value (an error exit, e.g. the -o output-conflict bug
            # which exits in seconds via HTCondor event 005, NOT 004). Both
            # are non-productive thrash. Everything else is a genuine 'ran'
            # attempt. end_reason is "was evicted" or "terminated".
            if dur_h is None:
                kind = "unknown"
            elif dur_h < IMMEDIATE_FAIL_H and (
                    a["end_reason"] == "was evicted"
                    or (a["end_reason"] == "terminated"
                        and a["retval"] not in (0, None))):
                kind = "immediate_fail"
            else:
                kind = "ran"
            rows.append({
                "leg": leg, "leg_type": leg_type,
                "campaign": "salvage" if is_salvage else "kartograf",
                "log": lp.name, "start": a["start"], "end": a["end"],
                "host": a["host"], "slot": a["slot"],
                "slot_type": classify_slot(a["slot"]),
                "device": a["device"], "capability": a["capability"],
                "end_reason": a["end_reason"], "retval": a["retval"],
                "duration_s": dur, "attempt_kind": kind,
                "leg_succeeded": leg_succeeded,
            })

    df = pd.DataFrame(rows)
    df["duration_h"] = df["duration_s"] / 3600.0
    df.to_csv(outdir / "compute_accounting_attempts.csv", index=False)
    print(f"Parsed {len(df)} execution attempts across "
          f"{df['leg'].nunique()} legs")

    # ---- Total burn ----
    total_h = df["duration_h"].sum()
    ran = df[df["attempt_kind"] == "ran"]
    imm = df[df["attempt_kind"] == "immediate_fail"]
    print("\n" + "=" * 64)
    print("TOTAL GPU BURN (all attempts)")
    print("=" * 64)
    print(f"  Total GPU-hours (all attempts):  {total_h:,.1f}")
    print(f"  Total attempts:                  {len(df):,}")
    print(f"  Total legs:                      {df['leg'].nunique():,}")
    print(f"  Mean attempts per leg:           "
          f"{len(df)/df['leg'].nunique():.1f}")

    # Attempt-level split: substantial run-attempts vs the seconds-long
    # -o output-conflict thrash. This is regime-independent: it just measures
    # GPU-hours in substantial vs trivial attempts, without claiming whether
    # substantial attempts accumulated (resume) or restarted (no-resume).
    print("\n  --- Attempt-level split (regime-independent) ---")
    print(f"  substantial attempts (>= 5 min):  {len(ran):,}  "
          f"({ran['duration_h'].sum():,.1f} GPU-h)")
    print(f"  trivial attempts (< 5 min evict): {len(imm):,}  "
          f"({imm['duration_h'].sum():,.1f} GPU-h)")
    print(f"  -> trivial -o-conflict thrash is {len(imm)/len(df)*100:.0f}% of "
          f"attempts but only {imm['duration_h'].sum()/total_h*100:.1f}% "
          f"of GPU-hours")
    print("\n  NOTE: Legs ran across preempting slots under two checkpointing")
    print("  regimes (--resume and no-resume versions of run_quickrun.sh), so")
    print("  we report total CONSUMED GPU-time, not productive-vs-redundant")
    print("  sampling. Per-leg 'time to complete' is intentionally not claimed.")

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
    # add genuine-compute (ran) GPU-hours per slot type
    ran_by_slot = (df[df["attempt_kind"] == "ran"].groupby("slot_type")
                   ["duration_h"].sum().rename("ran_gpu_hours"))
    slot = slot.merge(ran_by_slot, on="slot_type", how="left")
    slot["ran_gpu_hours"] = slot["ran_gpu_hours"].fillna(0)
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

    # ---- Per-ATTEMPT run-duration by GPU (regime-independent) ----
    # We report the duration distribution of substantial run-attempts
    # (>= 5 min), grouped by hardware. Each attempt is its own data point;
    # we do NOT aggregate to a per-leg "time to complete", which would
    # require resolving the resume/no-resume regime. This is an honest
    # "how long does a substantial simulation attempt take on each GPU".
    sub_att = df[df["attempt_kind"] == "ran"].copy()
    print("\n" + "=" * 64)
    print("SUBSTANTIAL RUN-ATTEMPT DURATION (per attempt, not per leg)")
    print("=" * 64)
    print(f"  Substantial run-attempts: {len(sub_att)}")
    if len(sub_att):
        print(f"  Mean: {sub_att['duration_h'].mean():.2f}h  "
              f"median: {sub_att['duration_h'].median():.2f}h  "
              f"min: {sub_att['duration_h'].min():.2f}h  "
              f"max: {sub_att['duration_h'].max():.2f}h")
        print("\n  By GPU device:")
        sd = sub_att.groupby("device").agg(
            n=("duration_h", "size"), mean_h=("duration_h", "mean"),
            median_h=("duration_h", "median"),
            min_h=("duration_h", "min"),
            max_h=("duration_h", "max")).reset_index()
        print(sd.sort_values("mean_h").to_string(index=False))
        print("\n  By capability:")
        sc = sub_att.groupby("capability").agg(
            n=("duration_h", "size"), mean_h=("duration_h", "mean"),
            median_h=("duration_h", "median")).reset_index()
        print(sc.sort_values("capability").to_string(index=False))
    success = sub_att  # used by the timeline/plots below

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

    # attempt-kind split: GPU-hours and attempt-counts (the -o churn story)
    ax = axes[1, 0]
    kinds = ["ran", "immediate_fail"]
    kh = [df[df["attempt_kind"] == k]["duration_h"].sum() for k in kinds]
    kn = [int((df["attempt_kind"] == k).sum()) for k in kinds]
    x = np.arange(len(kinds))
    ax2 = ax.twinx()
    ax.bar(x - 0.2, kh, width=0.4, color="#41ab5d", edgecolor="k",
           label="GPU-hours")
    ax2.bar(x + 0.2, kn, width=0.4, color="#fdae6b", edgecolor="k",
            label="attempt count")
    ax.set_xticks(x); ax.set_xticklabels(["ran\n(>=5min)",
                                          "immediate_fail\n(<5min evict)"])
    ax.set_ylabel("GPU-hours", color="#41ab5d")
    ax2.set_ylabel("attempt count", color="#e6550d")
    ax.set_title("-o conflict churn: hours vs attempt count")

    # successful timing by capability (mean bars)
    ax = axes[1, 1]
    if len(success):
        sc2 = success.groupby("capability")["duration_h"].mean()
        ax.bar(sc2.index.astype(str), sc2.values, color="#41b6c4",
               edgecolor="k")
        ax.set_xlabel("Capability")
        ax.set_ylabel("Mean run-attempt (h)")
        ax.set_title("Run-attempt duration by capability")

    # campaign split
    ax = axes[1, 2]
    ax.bar(camp["campaign"], camp["gpu_hours"], color=["#3182bd", "#e6550d"],
           edgecolor="k")
    ax.set_ylabel("GPU-hours"); ax.set_title("Burn: Kartograf vs salvage")

    plt.tight_layout()
    fig_path = outdir / "compute_accounting.png"
    plt.savefig(fig_path, dpi=130, bbox_inches="tight")
    print(f"\nWrote {fig_path}")

    # ---- Timeline: cumulative GPU-hours + concurrency over wall-clock ----
    # Concurrency uses RAN attempts only (genuine compute), so the peak
    # reflects real parallelism, not the -o conflict churn storms. We also
    # overlay all-attempts concurrency so the churn shows as the gap.
    tl_all = df.dropna(subset=["start", "end"]).copy()
    tl_all["start"] = pd.to_datetime(tl_all["start"])
    tl_all["end"] = pd.to_datetime(tl_all["end"])
    tl_ran = tl_all[tl_all["attempt_kind"] == "ran"].copy()

    def concurrency_curve(frame):
        events = []
        for _, r in frame.iterrows():
            events.append((r["start"], 1))
            events.append((r["end"], -1))
        ev = pd.DataFrame(events, columns=["t", "delta"]).sort_values("t")
        ev["concurrent"] = ev["delta"].cumsum()
        ev["dt_h"] = ev["t"].diff().dt.total_seconds().fillna(0) / 3600.0
        ev["cum_gpu_h"] = (ev["concurrent"].shift(1).fillna(0)
                           * ev["dt_h"]).cumsum()
        return ev

    if len(tl_ran):
        ev_ran = concurrency_curve(tl_ran)
        ev_all = concurrency_curve(tl_all)

        fig2, (axA, axB) = plt.subplots(2, 1, figsize=(14, 9), sharex=True)
        # all-attempts concurrency (churn-inflated) as faint background
        axA.fill_between(ev_all["t"], ev_all["concurrent"], alpha=0.18,
                         color="#969696", label="all attempts (incl. churn)")
        # ran-only concurrency (genuine parallelism)
        axA.plot(ev_ran["t"], ev_ran["concurrent"], color="#3182bd", lw=0.8)
        axA.fill_between(ev_ran["t"], ev_ran["concurrent"], alpha=0.35,
                         color="#6baed6", label="ran only (genuine compute)")
        peak_ran = int(ev_ran["concurrent"].max())
        peak_all = int(ev_all["concurrent"].max())
        axA.axhline(peak_ran, color="r", ls="--", lw=0.8,
                    label=f"peak ran = {peak_ran} GPUs")
        axA.set_ylabel("Concurrent GPUs in use")
        axA.set_title("Campaign concurrency over wall-clock time "
                      "(genuine compute vs churn)")
        axA.legend(fontsize=9)

        # cumulative GPU-hours: ran-only vs all
        axB.plot(ev_all["t"], ev_all["cum_gpu_h"], color="#969696", lw=1.2,
                 label=f"all attempts ({total_h:,.0f} GPU-h)")
        ran_total = tl_ran["duration_h"].sum()
        axB.plot(ev_ran["t"], ev_ran["cum_gpu_h"], color="#41ab5d", lw=1.5,
                 label=f"ran only ({ran_total:,.0f} GPU-h)")
        axB.set_ylabel("Cumulative GPU-hours")
        axB.set_xlabel("Wall-clock time")
        axB.set_title("Cumulative GPU burn: genuine compute vs total")
        axB.legend(fontsize=9)
        fig2.autofmt_xdate()
        plt.tight_layout()
        tl_path = outdir / "compute_timeline.png"
        plt.savefig(tl_path, dpi=130, bbox_inches="tight")
        print(f"Wrote {tl_path}")
        wall_span = (tl_all["end"].max() - tl_all["start"].min())
        print(f"\nWall-clock span:           {wall_span}")
        print(f"Peak concurrent (ran-only): {peak_ran} GPUs")
        print(f"Peak concurrent (all attempts, churn-inflated): {peak_all}")

    print(f"Wrote {outdir/'compute_accounting_attempts.csv'}")


if __name__ == "__main__":
    main()
