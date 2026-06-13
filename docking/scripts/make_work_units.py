#!/usr/bin/env python3
"""
docking/scripts/make_work_units.py  (v2 - one job per ligand)

Each work unit now corresponds to ONE cluster-representative ligand,
to be docked against ALL receptors (read from the shared receptors.zip
at job runtime, not enumerated here).

124 cluster representatives -> 124 work units -> 124 HTCondor jobs,
each docking 1 ligand x ~62 receptors sequentially
(~62-310 min/job at GNINA exhaustiveness=8).

Inputs:
  docking/ligands/cluster_representatives.csv
      columns: cluster_id, representative_test_compound, representative_smiles

Output:
  docking/work_units/lig_<cluster_id>.csv   (single-row CSV: cluster_id,
                                              ligand_name, ligand_smiles)
  docking/work_units/batch_list.txt         (list of work-unit filenames,
                                              for HTCondor queue-from-list)

Usage:
    python scripts/make_work_units.py
"""

import csv
from pathlib import Path

import pandas as pd

WORK_DIR = Path("docking/work_units")
WORK_DIR.mkdir(parents=True, exist_ok=True)


def main():
    ligs = pd.read_csv("docking/ligands/cluster_representatives.csv")
    print(f"Ligands (cluster representatives): {len(ligs)}")

    work_unit_names = []
    for _, row in ligs.iterrows():
        cid = int(row["cluster_id"])
        name = f"lig_{cid:04d}.csv"
        path = WORK_DIR / name
        with open(path, "w", newline="") as f:
            writer = csv.DictWriter(
                f, fieldnames=["cluster_id", "ligand_name", "ligand_smiles"])
            writer.writeheader()
            writer.writerow({
                "cluster_id": cid,
                "ligand_name": row["representative_test_compound"],
                "ligand_smiles": row["representative_smiles"],
            })
        work_unit_names.append(name)

    with open(WORK_DIR / "batch_list.txt", "w") as f:
        f.write("\n".join(work_unit_names) + "\n")

    print(f"\nWrote {len(work_unit_names)} work-unit files to {WORK_DIR}/")
    print(f"Wrote work-unit list to {WORK_DIR}/batch_list.txt")
    print("\nEach job docks 1 ligand against all receptors in "
          "receptors.zip (built by prep_receptors.py / "
          "zip_receptors step).")


if __name__ == "__main__":
    main()
