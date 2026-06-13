#!/usr/bin/env python3
"""
docking/scripts/make_work_units.py

Builds the full (cluster_representative x receptor) pair list and splits
it into batch CSVs of BATCH_SIZE rows each, for HTCondor submission.

124 cluster representatives x 64 receptors = 7,936 pairs.
At BATCH_SIZE=20, that's 397 batches, each ~20 dockings (~20-100 min
wall-clock per job at GNINA's default exhaustiveness=8).

Inputs:
  docking/ligands/cluster_representatives.csv
      columns: cluster_id, representative_test_compound, representative_smiles
  docking/receptors/boxes.csv
      columns: pdb_id, ligand_resname, center_x..size_z
      (produced by prep_receptors.py)

Output:
  docking/work_units/batch_NNNN.csv   (one per batch)
  docking/work_units/batch_list.txt   (list of batch filenames, for
                                        HTCondor queue-from-list)

Usage:
    python scripts/make_work_units.py [--batch-size 20]
"""

import argparse
import csv
from pathlib import Path

import pandas as pd

WORK_DIR = Path("docking/work_units")
WORK_DIR.mkdir(parents=True, exist_ok=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--batch-size", type=int, default=20,
                    help="Number of (ligand, receptor) dockings per "
                         "HTCondor job (default: 20)")
    args = ap.parse_args()

    ligs = pd.read_csv("docking/ligands/cluster_representatives.csv")
    boxes = pd.read_csv("docking/receptors/boxes.csv")

    print(f"Ligands (cluster representatives): {len(ligs)}")
    print(f"Receptors: {len(boxes)}")

    pdb_ids = boxes["pdb_id"].tolist()

    # Build the full cross-product
    rows = []
    for _, lrow in ligs.iterrows():
        for pdb_id in pdb_ids:
            rows.append({
                "cluster_id": lrow["cluster_id"],
                "ligand_name": lrow["representative_test_compound"],
                "ligand_smiles": lrow["representative_smiles"],
                "pdb_id": pdb_id,
            })

    print(f"Total (ligand, receptor) pairs: {len(rows)}")

    # Split into batches
    batch_size = args.batch_size
    n_batches = (len(rows) + batch_size - 1) // batch_size
    print(f"Batch size: {batch_size} -> {n_batches} batches")

    batch_names = []
    for i in range(n_batches):
        chunk = rows[i * batch_size:(i + 1) * batch_size]
        batch_name = f"batch_{i:04d}.csv"
        batch_path = WORK_DIR / batch_name
        with open(batch_path, "w", newline="") as f:
            writer = csv.DictWriter(
                f, fieldnames=["cluster_id", "ligand_name",
                               "ligand_smiles", "pdb_id"])
            writer.writeheader()
            writer.writerows(chunk)
        batch_names.append(batch_name)

    with open(WORK_DIR / "batch_list.txt", "w") as f:
        f.write("\n".join(batch_names) + "\n")

    print(f"\nWrote {n_batches} work-unit files to {WORK_DIR}/")
    print(f"Wrote batch list to {WORK_DIR}/batch_list.txt")
    print(f"\nEstimated wall-clock per batch at ~1-5 min/docking: "
          f"{batch_size*1}-{batch_size*5} minutes")


if __name__ == "__main__":
    main()
