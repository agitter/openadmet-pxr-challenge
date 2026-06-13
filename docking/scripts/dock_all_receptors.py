#!/usr/bin/env python3
"""
docking/scripts/dock_all_receptors.py

Dock one prepared ligand (PDBQT) against every receptor listed in
receptors/boxes.csv, using gnina (CNN rescoring enabled). Writes one
docked-pose SDF (gzipped) + log per receptor, plus a summary CSV.

Usage:
    python dock_all_receptors.py \
        --ligand-pdbqt ligand_3d/<cluster_id>.pdbqt \
        --boxes receptors/boxes.csv \
        --receptor-dir receptors \
        --cluster-id <cluster_id> \
        --ligand-name <ligand_name> \
        --outdir results/<cluster_id> \
        --summary results/<cluster_id>_summary.csv \
        --exhaustiveness 8 \
        --num-modes 5

boxes.csv columns expected:
    pdb_id,ligand_resname,center_x,center_y,center_z,size_x,size_y,size_z

For each row, looks for receptors/<pdb_id>_protein.pdbqt. Receptors
that can't be found are recorded as FAILED_NO_RECEPTOR in the summary
and skipped (does not abort the whole job).
"""

import argparse
import csv
import gzip
import shutil
import subprocess
import sys
from pathlib import Path


def dock_one(ligand_pdbqt, receptor_pdbqt, center, size, out_sdf, out_log,
              exhaustiveness, num_modes):
    cx, cy, cz = center
    sx, sy, sz = size
    cmd = [
        "gnina",
        "-r", str(receptor_pdbqt),
        "-l", str(ligand_pdbqt),
        "--center_x", str(cx), "--center_y", str(cy), "--center_z", str(cz),
        "--size_x", str(sx), "--size_y", str(sy), "--size_z", str(sz),
        "--exhaustiveness", str(exhaustiveness),
        "--num_modes", str(num_modes),
        "--cnn_scoring", "rescore",
        "-o", str(out_sdf),
        "--log", str(out_log),
    ]
    stdout_log = Path(str(out_log) + ".stdout")
    with open(stdout_log, "w") as f:
        result = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT)
    return result.returncode == 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ligand-pdbqt", required=True)
    ap.add_argument("--boxes", required=True,
                     help="Path to receptors/boxes.csv")
    ap.add_argument("--receptor-dir", required=True,
                     help="Directory containing <pdb_id>_protein.pdbqt files")
    ap.add_argument("--cluster-id", required=True)
    ap.add_argument("--ligand-name", required=True)
    ap.add_argument("--outdir", required=True,
                     help="Directory for per-receptor docked poses/logs")
    ap.add_argument("--summary", required=True,
                     help="Output summary CSV path")
    ap.add_argument("--exhaustiveness", type=int, default=8)
    ap.add_argument("--num-modes", type=int, default=5)
    args = ap.parse_args()

    ligand_pdbqt = Path(args.ligand_pdbqt)
    receptor_dir = Path(args.receptor_dir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if not ligand_pdbqt.exists():
        sys.exit(f"Ligand PDBQT not found: {ligand_pdbqt}")

    with open(args.boxes) as f:
        rows = list(csv.DictReader(f))

    summary_rows = []
    for row in rows:
        pdb_id = row["pdb_id"]
        receptor_pdbqt = receptor_dir / f"{pdb_id}_protein.pdbqt"
        tag = f"{args.cluster_id}__{pdb_id}"
        out_sdf = outdir / f"{tag}_docked.sdf"
        out_log = outdir / f"{tag}.log"

        if not receptor_pdbqt.exists():
            print(f"  [{tag}] receptor not found: {receptor_pdbqt}")
            summary_rows.append({
                "cluster_id": args.cluster_id,
                "ligand_name": args.ligand_name,
                "pdb_id": pdb_id,
                "status": "FAILED_NO_RECEPTOR",
            })
            continue

        center = (row["center_x"], row["center_y"], row["center_z"])
        size = (row["size_x"], row["size_y"], row["size_z"])

        print(f"--- Docking {args.ligand_name} (cluster {args.cluster_id}) "
              f"into {pdb_id} ---")
        ok = dock_one(ligand_pdbqt, receptor_pdbqt, center, size,
                       out_sdf, out_log, args.exhaustiveness, args.num_modes)

        if ok and out_sdf.exists():
            gz_path = Path(str(out_sdf) + ".gz")
            with open(out_sdf, "rb") as f_in, gzip.open(gz_path, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
            out_sdf.unlink()
            print(f"  [{tag}] done")
            summary_rows.append({
                "cluster_id": args.cluster_id,
                "ligand_name": args.ligand_name,
                "pdb_id": pdb_id,
                "status": "success",
            })
        else:
            print(f"  [{tag}] gnina FAILED, see {out_log}.stdout")
            summary_rows.append({
                "cluster_id": args.cluster_id,
                "ligand_name": args.ligand_name,
                "pdb_id": pdb_id,
                "status": "FAILED",
            })

    with open(args.summary, "w", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=["cluster_id", "ligand_name", "pdb_id", "status"])
        writer.writeheader()
        writer.writerows(summary_rows)

    n_ok = sum(1 for r in summary_rows if r["status"] == "success")
    print(f"\n{n_ok}/{len(summary_rows)} receptors docked successfully")
    print(f"Summary written to {args.summary}")


if __name__ == "__main__":
    main()
