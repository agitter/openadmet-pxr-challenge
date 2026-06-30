#!/usr/bin/env python3
"""
diagnose_slots.py

Diagnose the slot-type classification by examining the actual SlotName
prefixes on each host. Reveals whether researcher-owned CHTC nodes are
falling into the 'shared' bucket when they should be backfill (or
whatever the correct category is).

Run from repo root after 22_compute_accounting.py has written
openfe/compute_accounting_attempts.csv.

Usage:
    python diagnose_slots.py
"""

import pandas as pd

df = pd.read_csv("openfe/compute_accounting_attempts.csv")

# Strip trailing digits/underscores from the slot prefix for grouping
df["slot_prefix"] = (df["slot"].str.split("@").str[0]
                     .str.replace(r"[0-9_]+$", "", regex=True))

print("=" * 64)
print("1. SLOT PREFIX vs ASSIGNED slot_type")
print("=" * 64)
print(df.groupby(["slot_prefix", "slot_type"]).size()
      .sort_values(ascending=False).head(40).to_string())

print("\n" + "=" * 64)
print("2. HOSTS CURRENTLY CLASSIFIED AS 'shared' (top 25 by attempts)")
print("=" * 64)
shared = df[df["slot_type"] == "shared"]
print(shared["host"].value_counts().head(25).to_string())

print("\n" + "=" * 64)
print("3. RAW SLOT NAMES IN 'shared' BUCKET (distinct samples)")
print("=" * 64)
print(shared["slot"].dropna().drop_duplicates().head(25).to_string())

print("\n" + "=" * 64)
print("4. SLOT PREFIX BREAKDOWN PER REPRESENTATIVE HOST")
print("=" * 64)
print("(reveals whether a host uses backfill* or slotN_M* prefixes)")
check_hosts = [
    "jcaicedogpu0002.chtc.wisc.edu", "vetsigian0001.chtc.wisc.edu",
    "ahlquist0000.chtc.wisc.edu", "blengerichgpu4000.chtc.wisc.edu",
    "mkhodakgpu4000.chtc.wisc.edu", "amuraligpu4000.chtc.wisc.edu",
    "gpulab2003.chtc.wisc.edu", "gpulab2005.chtc.wisc.edu",
    "gpu4002.chtc.wisc.edu", "gpu4006.chtc.wisc.edu",
    "gpu2008.chtc.wisc.edu", "gpu2011.chtc.wisc.edu",
    "gitter0000.chtc.wisc.edu", "gpulab2001.chtc.wisc.edu",
]
for h in check_hosts:
    sub = df[df["host"] == h]
    if len(sub):
        prefixes = dict(sub["slot"].str.split("@").str[0]
                        .str.replace(r"[0-9_]+$", "", regex=True)
                        .value_counts())
        types = dict(sub["slot_type"].value_counts())
        mean_dur = sub["duration_h"].mean()
        print(f"\n{h}")
        print(f"  prefixes:   {prefixes}")
        print(f"  slot_types: {types}")
        print(f"  mean attempt duration: {mean_dur:.3f}h  "
              f"(n attempts={len(sub)})")

print("\n" + "=" * 64)
print("5. ALL DISTINCT (slot_prefix, host-domain) COMBINATIONS")
print("=" * 64)
df["host_domain"] = df["host"].str.extract(r"\.([\w]+\.[\w]+)$")[0].fillna(
    "short-name")
combo = (df.groupby(["slot_prefix", "host_domain"])
         .agg(attempts=("duration_h", "size"),
              mean_dur_h=("duration_h", "mean"))
         .reset_index().sort_values("attempts", ascending=False))
print(combo.head(30).to_string(index=False))
