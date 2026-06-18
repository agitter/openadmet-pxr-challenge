#!/usr/bin/env python3
"""
openfe/scripts/04_prepare_planning_inputs.py

For each of the 124 clusters, assemble a single multi-molecule SDF
containing all ligands in that cluster's network (test compounds +
training anchors), ready for `openfe plan-rbfe-network -M <sdf>`.

The per-cluster SDF is built by concatenating the individual
<ligand_id>_ligand.sdf files from openfe/rbfe_inputs/<cluster_id>/.

Also writes:
  openfe/plan_work_units/cluster_list.txt
      one line per cluster_id, used as the HTCondor queue source

Usage:
    pip install pandas rdkit
    python openfe/scripts/04_prepare_planning_inputs.py \
        --rbfe-inputs openfe/rbfe_inputs \
        --outdir openfe/plan_inputs \
        --work-units openfe/plan_work_units
"""

import argparse
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rbfe-inputs", default="openfe/rbfe_inputs",
                    help="Directory containing per-cluster subdirectories "
                         "with <ligand_id>_ligand.sdf files")
    ap.add_argument("--outdir", default="openfe/plan_inputs",
                    help="Output directory for per-cluster SDFs and "
                         "protein PDB symlinks")
    ap.add_argument("--work-units", default="openfe/plan_work_units",
                    help="Output directory for HTCondor work-unit files")
    args = ap.parse_args()

    rbfe_inputs = Path(args.rbfe_inputs)
    outdir = Path(args.outdir)
    work_units = Path(args.work_units)
    outdir.mkdir(parents=True, exist_ok=True)
    work_units.mkdir(parents=True, exist_ok=True)

    # Find all cluster subdirectories
    cluster_dirs = sorted(
        [d for d in rbfe_inputs.iterdir() if d.is_dir()
         and d.name != "extraction_summary_extended.csv"],
        key=lambda d: int(d.name)
    )
    print(f"Found {len(cluster_dirs)} cluster directories")

    summary_rows = []
    cluster_ids = []

    for cluster_dir in cluster_dirs:
        cluster_id = cluster_dir.name

        # Find all ligand SDF files in this cluster directory
        ligand_sdfs = sorted(cluster_dir.glob("*_ligand.sdf"))
        if not ligand_sdfs:
            print(f"  WARNING: no ligand SDFs found in {cluster_dir}")
            continue

        # Find the receptor PDB (all should be the same in a cluster;
        # take the first one - they're copies of the same prepared PDB)
        receptor_pdbs = sorted(cluster_dir.glob("*_receptor.pdb"))
        if not receptor_pdbs:
            print(f"  WARNING: no receptor PDB found in {cluster_dir}")
            continue
        receptor_pdb = receptor_pdbs[0]

        # Read and concatenate all ligand mols into one SDF
        cluster_out_dir = outdir / cluster_id
        cluster_out_dir.mkdir(exist_ok=True)
        combined_sdf = cluster_out_dir / "ligands.sdf"
        combined_receptor = cluster_out_dir / "receptor.pdb"

        mols = []
        failed = []
        for sdf in ligand_sdfs:
            suppl = Chem.SDMolSupplier(str(sdf), removeHs=False,
                                        sanitize=False)
            for mol in suppl:
                if mol is None:
                    failed.append(str(sdf))
                    continue
                try:
                    Chem.SanitizeMol(
                        mol,
                        sanitizeOps=Chem.SANITIZE_ALL
                                     ^ Chem.SANITIZE_KEKULIZE)
                except Exception:
                    pass
                mols.append(mol)

        if not mols:
            print(f"  {cluster_id}: ERROR - no valid mols read from "
                  f"{len(ligand_sdfs)} SDFs")
            continue

        # Write combined SDF
        writer = Chem.SDWriter(str(combined_sdf))
        writer.SetKekulize(False)
        for mol in mols:
            writer.write(mol)
        writer.close()

        # Copy receptor PDB (rename to receptor.pdb for simplicity)
        import shutil
        shutil.copy(receptor_pdb, combined_receptor)

        n_test = sum(1 for s in ligand_sdfs if s.stem.startswith("T"))
        n_anchor = sum(1 for s in ligand_sdfs if s.stem.startswith("A"))

        summary_rows.append({
            "cluster_id": cluster_id,
            "n_ligands": len(mols),
            "n_test": n_test,
            "n_anchors": n_anchor,
            "n_failed_mols": len(failed),
            "ligands_sdf": str(combined_sdf),
            "receptor_pdb": str(combined_receptor),
        })
        cluster_ids.append(cluster_id)

        if len(cluster_ids) % 20 == 0 or len(cluster_ids) == 1:
            print(f"  Prepared cluster {cluster_id} "
                  f"({len(cluster_ids)}/{len(cluster_dirs)}): "
                  f"{len(mols)} ligands ({n_test} test + {n_anchor} anchor)")

    # Write cluster list for HTCondor queue
    cluster_list = work_units / "cluster_list.txt"
    with open(cluster_list, "w") as f:
        f.write("\n".join(cluster_ids) + "\n")

    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(outdir / "planning_inputs_summary.csv", index=False)

    print(f"\nPrepared {len(cluster_ids)} / {len(cluster_dirs)} clusters")
    print(f"Wrote cluster list -> {cluster_list}")
    print(f"Wrote summary -> {outdir/'planning_inputs_summary.csv'}")
    print(f"\nLigands per cluster:")
    print(summary_df["n_ligands"].describe().to_string())


if __name__ == "__main__":
    main()
