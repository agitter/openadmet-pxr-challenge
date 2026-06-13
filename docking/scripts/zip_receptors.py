#!/usr/bin/env python3
"""
docking/scripts/zip_receptors.py

Bundles all receptor PDBQT files + boxes.csv into a single
docking/receptors.zip, so HTCondor transfers ONE file per job instead
of ~62 small files.

Run after prep_receptors.py.

Usage:
    python scripts/zip_receptors.py
"""

import zipfile
from pathlib import Path

RECEPTOR_DIR = Path("docking/receptors")
OUT_ZIP = Path("docking/receptors.zip")


def main():
    pdbqt_files = sorted(RECEPTOR_DIR.glob("*_protein.pdbqt"))
    boxes_csv = RECEPTOR_DIR / "boxes.csv"

    if not pdbqt_files:
        raise SystemExit(f"No *_protein.pdbqt files found in {RECEPTOR_DIR} "
                          f"- run prep_receptors.py first.")
    if not boxes_csv.exists():
        raise SystemExit(f"{boxes_csv} not found - run prep_receptors.py first.")

    with zipfile.ZipFile(OUT_ZIP, "w", zipfile.ZIP_DEFLATED) as zf:
        for f in pdbqt_files:
            zf.write(f, arcname=f.name)
        zf.write(boxes_csv, arcname=boxes_csv.name)

    print(f"Wrote {OUT_ZIP} containing {len(pdbqt_files)} receptor PDBQTs "
          f"+ boxes.csv")
    print(f"Size: {OUT_ZIP.stat().st_size / 1e6:.1f} MB")


if __name__ == "__main__":
    main()
