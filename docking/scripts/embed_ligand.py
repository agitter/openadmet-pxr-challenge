#!/usr/bin/env python3
"""
docking/scripts/embed_ligand.py

Convert a single SMILES string to a 3D-embedded SDF file using RDKit
(ETKDGv3 embedding + MMFF94 optimization, falling back to UFF if MMFF
parameters are unavailable for the given molecule).

Usage:
    python embed_ligand.py "<SMILES>" <output.sdf>

Exits non-zero (with a message to stderr) on parse or embedding
failure, so callers (run_batch.sh) can detect and report errors per
ligand.
"""

import sys

from rdkit import Chem
from rdkit.Chem import AllChem


def main():
    if len(sys.argv) != 3:
        sys.exit("Usage: python embed_ligand.py <SMILES> <output.sdf>")

    smiles, out_path = sys.argv[1], sys.argv[2]

    # Strip CXSMILES extended notation if present (e.g. "... |&1:7,9|" suffix
    # encoding enhanced stereochemistry with "either" stereocenters). RDKit's
    # MolFromSmiles does not handle CXSMILES syntax; stripping the | extension
    # loses the either-stereochemistry annotation but preserves the core
    # structure and defined stereocenters, which is sufficient for 3D embedding.
    if " |" in smiles:
        smiles = smiles[:smiles.index(" |")]

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        sys.exit(f"Failed to parse SMILES: {smiles}")

    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42

    if AllChem.EmbedMolecule(mol, params) < 0:
        sys.exit(f"3D embedding failed for SMILES: {smiles}")

    try:
        AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
    except Exception:
        AllChem.UFFOptimizeMolecule(mol, maxIters=500)

    writer = Chem.SDWriter(out_path)
    writer.write(mol)
    writer.close()


if __name__ == "__main__":
    main()
