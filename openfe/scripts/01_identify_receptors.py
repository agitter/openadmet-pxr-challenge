#!/usr/bin/env python3
"""
openfe/scripts/01_identify_receptors.py

Identify the unique PXR receptor structures needed for the RBFE
campaign, based on rbfe_template_selection.csv (the per-cluster
receptor/pose selection from the ensemble docking).

Reports:
  - How many unique receptors are needed (from the 62-structure ensemble)
  - How many clusters use each receptor (usage count)
  - Which receptor PDB files need to be prepared with PDBFixer

Writes:
  openfe/receptors/receptor_usage.csv
      columns: pdb_id, n_clusters, cluster_ids
      (one row per unique receptor, sorted by usage count descending)

Usage:
    python openfe/scripts/01_identify_receptors.py \
        --selection docking/docking_analysis/rbfe_template_selection.csv \
        --outdir openfe/receptors
"""

import argparse
from pathlib import Path

import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--selection", required=True,
                    help="rbfe_template_selection.csv")
    ap.add_argument("--outdir", default="openfe/receptors")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    sel = pd.read_csv(args.selection)
    print(f"Loaded {len(sel)} cluster selections")
    print(f"Unique receptors needed: {sel['selected_pdb_id'].nunique()} "
          f"of 62 available")

    usage = (sel.groupby("selected_pdb_id")
               .agg(
                   n_clusters=("cluster_id", "count"),
                   cluster_ids=("cluster_id", lambda x: ",".join(str(v) for v in sorted(x))),
               )
               .reset_index()
               .rename(columns={"selected_pdb_id": "pdb_id"})
               .sort_values("n_clusters", ascending=False)
               .reset_index(drop=True))

    print("\nReceptor usage (sorted by cluster count):")
    print(usage[["pdb_id", "n_clusters"]].to_string(index=False))

    out_path = outdir / "receptor_usage.csv"
    usage.to_csv(out_path, index=False)
    print(f"\nWrote {out_path}")
    print("\nNext step: run PDBFixer on each of these receptors.")
    print("Input PDBs are at:")
    print("  external/pxr_xtal_re-refinement/pxr_rerefined_structures/"
          "<pdb_id>/<pdb_id>.pdb")


if __name__ == "__main__":
    main()
