#!/usr/bin/env python3
"""
docking/scripts/make_work_units_extended.py

Generate HTCondor work units for the extended docking campaign (all
602 ligands: 513 test compounds + 89 unique training anchors, each
against all ~62 receptors).

Produces the SAME single-row CSV schema as the original
make_work_units.py (cluster_id, ligand_name, ligand_smiles), so
run_batch.sh / embed_ligand.py / dock_all_receptors.py / receptors.zip
are reused UNCHANGED - only the "cluster_id" column now holds a string
ligand_id (e.g. "T042" or "A0063") instead of an integer cluster
label. These ligand_id values are filesystem-safe and do not collide
with the original campaign's purely-numeric cluster_id directories
(e.g. "21", "401") already present under results/.

Inputs:
  docking/all_rbfe_ligands.csv (from build_extended_ligand_list.py)
      columns: ligand_id, ligand_name, smiles, source

Output:
  docking/work_units/lig_<ligand_id>.csv   (602 files)
  docking/work_units/batch_list_extended.txt

Usage:
    python make_work_units_extended.py
"""

import csv
from pathlib import Path

import pandas as pd

WORK_DIR = Path("docking/work_units")
WORK_DIR.mkdir(parents=True, exist_ok=True)


def main():
    ligands = pd.read_csv("docking/all_rbfe_ligands.csv")
    print(f"Ligands: {len(ligands)}")

    work_unit_names = []
    for _, row in ligands.iterrows():
        ligand_id = row["ligand_id"]
        name = f"lig_{ligand_id}.csv"
        path = WORK_DIR / name
        with open(path, "w", newline="") as f:
            writer = csv.DictWriter(
                f, fieldnames=["cluster_id", "ligand_name", "ligand_smiles"])
            writer.writeheader()
            writer.writerow({
                "cluster_id": ligand_id,
                "ligand_name": row["ligand_name"],
                "ligand_smiles": row["smiles"],
            })
        work_unit_names.append(name)

    with open(WORK_DIR / "batch_list_extended.txt", "w") as f:
        f.write("\n".join(work_unit_names) + "\n")

    print(f"\nWrote {len(work_unit_names)} work-unit files to {WORK_DIR}/")
    print(f"Wrote work-unit list to {WORK_DIR}/batch_list_extended.txt")
    print("\nEach job docks 1 ligand against all receptors in "
          "receptors.zip (already built from the original campaign - "
          "no need to rebuild).")
    print("\nSubmit with the EXISTING submit_docking.sub, overriding the "
          "queue source:")
    print('  condor_submit submit_docking.sub '
          '-append "queue work_unit from work_units/batch_list_extended.txt"')


if __name__ == "__main__":
    main()
