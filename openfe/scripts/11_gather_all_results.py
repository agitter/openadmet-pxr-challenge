#!/usr/bin/env python3
"""
openfe/scripts/11_gather_all_results.py

Gather ALL RBFE leg results into a single authoritative per-edge table,
combining the original Kartograf run and the LOMAP salvage run, with
provenance tracking.

Result locations searched per leg (in priority order):
  1. <job>/salvage/result.json        - LOMAP salvage (mapper=lomap)
  2. <job>/result.json                - rerun without checkpoint (kartograf)
  3. <job>/quickrun_output/result.json - original checkpoint run (kartograf)

For each edge (complex + solvent pair), computes ddG = dG_complex -
dG_solvent and records which mapper produced each leg. A salvaged leg is
only used if the original Kartograf leg failed (we never override a good
Kartograf result with a salvage).

Provenance tiers:
  - tier1_kartograf: both legs from Kartograf
  - tier2_lomap: at least one leg from LOMAP salvage
  - incomplete: at least one leg has no valid result anywhere

Outputs:
  - openfe/all_edge_results.csv: per-edge ddG, mapper tier, leg sources
  - openfe/docking_fallback_compounds.csv: compounds with no usable edge

Usage:
    python openfe/scripts/11_gather_all_results.py \
        --production-dir openfe/production \
        --outdir openfe
"""

import argparse
import json
import re
from collections import defaultdict
from pathlib import Path

import pandas as pd


def load_estimate(result_path):
    """Return estimate magnitude (float) or None."""
    if not result_path.exists() or result_path.stat().st_size <= 100:
        return None
    try:
        d = json.loads(result_path.read_text())
        est = d.get("estimate")
        if est and est.get("magnitude") is not None:
            return est["magnitude"]
    except Exception:
        pass
    return None


def get_leg_result(job_dir):
    """Find a valid result for a leg, checking all locations.
    Returns (dG, mapper) or (None, None)."""
    # Salvage (LOMAP) first - but only matters if Kartograf failed,
    # which is handled by the caller. Here we just report what exists.
    salvage = load_estimate(job_dir / "salvage" / "result.json")
    top = load_estimate(job_dir / "result.json")
    inner = load_estimate(job_dir / "quickrun_output" / "result.json")

    # Prefer Kartograf (top or inner) if valid; fall back to salvage.
    if top is not None:
        return top, "kartograf"
    if inner is not None:
        return inner, "kartograf"
    if salvage is not None:
        return salvage, "lomap"
    return None, None


def edge_and_leg(transform_name):
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

    # edge -> {complex: (dG, mapper), solvent: (dG, mapper)}
    edges = defaultdict(dict)

    with open(transform_list) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            cluster_id, transform_name = [x.strip() for x in line.split(",", 1)]
            job_dir = prod_dir / cluster_id / transform_name
            edge, leg = edge_and_leg(transform_name)
            if leg is None:
                continue
            dG, mapper = get_leg_result(job_dir)
            edges[(cluster_id, edge)][leg] = (dG, mapper)

    rows = []
    for (cluster_id, edge), legs in edges.items():
        c_dG, c_mapper = legs.get("complex", (None, None))
        s_dG, s_mapper = legs.get("solvent", (None, None))

        both_done = c_dG is not None and s_dG is not None
        ddg = (c_dG - s_dG) if both_done else None

        # Determine tier
        mappers = {m for m in [c_mapper, s_mapper] if m}
        if not both_done:
            tier = "incomplete"
        elif mappers == {"kartograf"}:
            tier = "tier1_kartograf"
        elif "lomap" in mappers:
            tier = "tier2_lomap"
        else:
            tier = "unknown"

        rows.append({
            "cluster_id": cluster_id,
            "edge": edge,
            "both_done": both_done,
            "ddg": ddg,
            "tier": tier,
            "complex_dG": c_dG,
            "complex_mapper": c_mapper,
            "solvent_dG": s_dG,
            "solvent_mapper": s_mapper,
        })

    df = pd.DataFrame(rows)
    out_path = outdir / "all_edge_results.csv"
    df.to_csv(out_path, index=False)

    print("=" * 60)
    print("COMBINED RESULTS (Kartograf + LOMAP salvage)")
    print("=" * 60)
    print(f"Total edges: {len(df)}")
    print(f"Both legs done: {df['both_done'].sum()}")
    print(f"Incomplete: {(~df['both_done']).sum()}")
    print()
    print("By tier:")
    print(df["tier"].value_counts().to_string())

    done = df[df["both_done"]]
    if len(done):
        print(f"\nddG distribution (all usable edges):")
        ddgs = done["ddg"]
        print(f"  mean={ddgs.mean():.2f}  median={ddgs.median():.2f}  "
              f"min={ddgs.min():.2f}  max={ddgs.max():.2f}")
        n_salvage = (done["tier"] == "tier2_lomap").sum()
        print(f"\nEdges recovered by LOMAP salvage: {n_salvage}")
        print(f"Original Kartograf edges: {(done['tier']=='tier1_kartograf').sum()}")

    print(f"\nWrote {out_path}")


if __name__ == "__main__":
    main()
