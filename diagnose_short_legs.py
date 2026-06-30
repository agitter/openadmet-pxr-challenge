#!/usr/bin/env python3
"""
diagnose_short_legs.py

Identify why some "successful legs" have implausibly short productive runs
(e.g. 7 seconds). Hypothesis: these legs completed via checkpoint-resume -
most work happened in an earlier (evicted) attempt, and the final attempt
that produced result.json just resumed and finished in minutes. If so, the
per-leg metric should be the SUM of ran attempts, not the longest single one.

Run from repo root after 22_compute_accounting.py wrote the attempts CSV.

Usage:
    python diagnose_short_legs.py
"""

import pandas as pd

df = pd.read_csv("openfe/compute_accounting_attempts.csv")

# Reproduce the current success selection: legs with result.json,
# ran attempts, longest ran attempt per leg.
succ = df[df["leg_succeeded"] & (df["attempt_kind"] == "ran")].copy()
idx = succ.groupby("leg")["duration_h"].idxmax()
chosen = succ.loc[idx]

tiny = chosen[chosen["duration_h"] < 0.1].sort_values("duration_h")
print("=" * 64)
print(f"Legs whose LONGEST ran attempt is < 0.1h (6 min): {len(tiny)}")
print("=" * 64)
print(tiny[["leg", "duration_h", "end_reason", "retval", "host",
            "device"]].head(20).to_string(index=False))

# For the shortest few, show ALL attempts to understand the structure
print("\n" + "=" * 64)
print("FULL ATTEMPT STRUCTURE FOR 3 SHORTEST 'SUCCESSFUL' LEGS")
print("=" * 64)
for bad_leg in tiny["leg"].head(3):
    allatt = df[df["leg"] == bad_leg].sort_values("duration_h",
                                                  ascending=False)
    ran = allatt[allatt["attempt_kind"] == "ran"]
    print(f"\n--- leg: {bad_leg} ---")
    print(allatt[["log", "duration_h", "end_reason", "retval",
                  "attempt_kind", "host"]].head(15).to_string(index=False))
    print(f"  total attempts: {len(allatt)}  | ran attempts: {len(ran)}")
    print(f"  longest single ran attempt: {ran['duration_h'].max():.4f}h")
    print(f"  SUM of ran attempts:        {ran['duration_h'].sum():.4f}h")

# Compare the two candidate per-leg metrics across ALL successful legs
print("\n" + "=" * 64)
print("LONGEST-ATTEMPT vs SUM-OF-RAN-ATTEMPTS (all successful legs)")
print("=" * 64)
per_leg = succ.groupby("leg")["duration_h"].agg(
    longest="max", sum_ran="sum", n_ran="size").reset_index()
print(f"  legs: {len(per_leg)}")
print(f"  longest-attempt metric: mean={per_leg['longest'].mean():.2f}h  "
      f"median={per_leg['longest'].median():.2f}h  "
      f"min={per_leg['longest'].min():.4f}h")
print(f"  sum-of-ran metric:      mean={per_leg['sum_ran'].mean():.2f}h  "
      f"median={per_leg['sum_ran'].median():.2f}h  "
      f"min={per_leg['sum_ran'].min():.4f}h")
print(f"\n  legs with >1 ran attempt (fragmented by preemption): "
      f"{(per_leg['n_ran'] > 1).sum()}")
print(f"  legs with exactly 1 ran attempt: {(per_leg['n_ran'] == 1).sum()}")

# How many tiny-longest legs become reasonable under sum?
merged = per_leg[per_leg["longest"] < 0.1]
print(f"\n  Of the {len(merged)} legs with longest<0.1h:")
print(f"    their SUM-of-ran: mean={merged['sum_ran'].mean():.2f}h  "
      f"min={merged['sum_ran'].min():.4f}h  "
      f"max={merged['sum_ran'].max():.2f}h")
print(f"    still <0.1h under sum metric: "
      f"{(merged['sum_ran'] < 0.1).sum()}")
