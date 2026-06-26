#!/usr/bin/env python3
"""
openfe/scripts/10_salvage_plan.py

Re-plan FAILED edges using the LOMAP atom mapper (element_change=False)
instead of Kartograf. Salvages only edges that failed under Kartograf,
preserving all original results for reproducibility and compute accounting.

Output layout: each salvage attempt goes into a NEW 'salvage/' subdirectory
inside the existing job directory, so the original (failed) Kartograf
result.json, logs, and quickrun_output are never touched:

  openfe/production/<cluster_id>/<transform_name>/
    result.json              <- original Kartograf result (KEPT, may be null)
    quickrun_output/         <- original Kartograf data (KEPT)
    logs/                    <- original Kartograf logs (KEPT)
    <transform_name>.json    <- original Kartograf transformation (KEPT)
    salvage/                 <- NEW salvage attempt, isolated
      <transform_name>.json  <- LOMAP transformation JSON
      run_quickrun.sh        <- copied executable
      (result.json written here by the production run)

Provenance: each salvage transformation JSON is tagged with mapper info
(mapper=lomap, element_change=False, max3d, n_mapped_atoms) so the final
pEC50 analysis can treat salvaged edges as a separate confidence tier.

Must run INSIDE the openfe container (needs openfe Python API).

Usage (inside container, on AP) - process in batches of 50:
    # Batch 1 (edges 0-49, truncates output files):
    python openfe/scripts/10_salvage_plan.py \
        --edge-report openfe/final_edge_report.csv \
        --production-dir openfe/production \
        --rbfe-inputs openfe/rbfe_inputs \
        --script openfe/scripts/run_quickrun.sh \
        --batch-start 0 --batch-size 50

    # Batch 2 (edges 50-99, appends), then --batch-start 100, 150, 200...
    # The script prints the next --batch-start to use.

All edges with both_done=False are salvaged (every incomplete edge).
After all batches, re-run connectivity analysis on the combined
Kartograf + LOMAP results to see total compound coverage.
"""

import argparse
import json
from pathlib import Path

import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--edge-report", default="openfe/final_edge_report.csv",
                    help="final_edge_report.csv; all edges with "
                         "both_done=False are salvaged")
    ap.add_argument("--production-dir", default="openfe/production")
    ap.add_argument("--rbfe-inputs", default="openfe/rbfe_inputs")
    ap.add_argument("--script", default="openfe/scripts/run_quickrun.sh")
    ap.add_argument("--min-mapped-frac", type=float, default=0.25,
                    help="Skip edges where the mapping covers less than "
                         "this fraction of the smaller ligand's heavy "
                         "atoms (default 0.25). Low-coverage mappings give "
                         "unreliable ddG and often NaN.")
    ap.add_argument("--min-mapped-atoms", type=int, default=4,
                    help="Also require at least this many mapped atoms "
                         "(default 4).")
    ap.add_argument("--outlist", default="openfe/production/salvage_transform_list.txt",
                    help="Output queue list for the salvage production run")
    ap.add_argument("--batch-start", type=int, default=0,
                    help="Index of first failed edge to process (0-based). "
                         "For batching: 0, 50, 100, ... (default 0)")
    ap.add_argument("--batch-size", type=int, default=50,
                    help="Max edges to process this invocation (default 50). "
                         "Keeps each AP run short.")
    args = ap.parse_args()

    import shutil
    import openfe
    from openfe import (SmallMoleculeComponent, ProteinComponent,
                        SolventComponent, ChemicalSystem)
    from openfe.setup import LomapAtomMapper
    from openfe.protocols.openmm_rfe import RelativeHybridTopologyProtocol
    from gufe import Transformation
    from rdkit import Chem
    from rdkit import RDLogger
    RDLogger.DisableLog("rdApp.*")

    prod_dir = Path(args.production_dir)
    rbfe_inputs = Path(args.rbfe_inputs)
    script_src = Path(args.script)

    salvage = pd.read_csv(args.edge_report)
    salvage = salvage[~salvage["both_done"]].copy().reset_index(drop=True)
    total_failed = len(salvage)

    # Slice to the requested batch
    batch = salvage.iloc[args.batch_start:args.batch_start + args.batch_size]
    print(f"Total failed edges: {total_failed}")
    print(f"Processing batch: edges {args.batch_start} to "
          f"{args.batch_start + len(batch) - 1} ({len(batch)} edges)")

    def find_ligand_file(cluster_id, mol_name):
        cluster_dir = rbfe_inputs / str(cluster_id)
        for sdf in cluster_dir.glob("*_ligand.sdf"):
            suppl = Chem.SDMolSupplier(str(sdf), removeHs=False)
            for mol in suppl:
                if mol and mol.GetProp("_Name") == mol_name:
                    return sdf
        return None

    # LOMAP mapper: element changes forbidden, slightly looser 3D cutoff
    lomap_mapper = LomapAtomMapper(
        threed=True, max3d=1.5, element_change=False)

    settings = RelativeHybridTopologyProtocol.default_settings()
    settings.protocol_repeats = 1
    protocol = RelativeHybridTopologyProtocol(settings)
    solvent = SolventComponent()

    transform_entries = []
    n_ok = n_no_mapping = n_error = n_low_coverage = 0
    skipped_edges = []

    for _, row in batch.iterrows():
        cluster_id = row["cluster_id"]
        edge = row["edge"]
        body = edge[len("rbfe_"):] if edge.startswith("rbfe_") else edge
        parts = body.split("_")
        if len(parts) != 2:
            print(f"  SKIP {edge}: cannot parse ligand names")
            n_error += 1
            continue
        nameA, nameB = parts

        try:
            sdfA = find_ligand_file(cluster_id, nameA)
            sdfB = find_ligand_file(cluster_id, nameB)
            if sdfA is None or sdfB is None:
                print(f"  SKIP {cluster_id}/{edge}: ligand SDF not found")
                n_error += 1
                continue

            molA = SmallMoleculeComponent.from_sdf_file(str(sdfA))
            molB = SmallMoleculeComponent.from_sdf_file(str(sdfB))
            receptor_pdb = next((rbfe_inputs / str(cluster_id)).glob(
                "*_receptor.pdb"))
            protein = ProteinComponent.from_pdb_file(str(receptor_pdb))

            mappings = list(lomap_mapper.suggest_mappings(molA, molB))
            if not mappings:
                print(f"  NO MAPPING {cluster_id}/{edge}: LOMAP also failed")
                n_no_mapping += 1
                continue
            mapping = mappings[0]
            n_mapped = len(mapping.componentA_to_componentB)

            # Quality check: reject low-coverage mappings that would give
            # unreliable ddG / likely NaN. Use the smaller ligand's heavy
            # atom count as the denominator.
            nheavy_A = molA.to_rdkit().GetNumHeavyAtoms()
            nheavy_B = molB.to_rdkit().GetNumHeavyAtoms()
            smaller = min(nheavy_A, nheavy_B)
            frac = n_mapped / smaller if smaller > 0 else 0.0

            if n_mapped < args.min_mapped_atoms or frac < args.min_mapped_frac:
                print(f"  SKIP {cluster_id}/{edge}: low coverage "
                      f"({n_mapped} atoms = {frac:.0%} of smaller ligand, "
                      f"below thresholds) - unsalvageable, will use docking")
                n_low_coverage += 1
                skipped_edges.append({
                    "cluster_id": cluster_id, "edge": edge,
                    "n_mapped": n_mapped, "frac": round(frac, 3),
                    "reason": "low_coverage",
                })
                continue

            sysA_c = ChemicalSystem(
                {"ligand": molA, "protein": protein, "solvent": solvent})
            sysB_c = ChemicalSystem(
                {"ligand": molB, "protein": protein, "solvent": solvent})
            sysA_s = ChemicalSystem({"ligand": molA, "solvent": solvent})
            sysB_s = ChemicalSystem({"ligand": molB, "solvent": solvent})

            for stateA, stateB, leg in [
                (sysA_c, sysB_c, "complex"),
                (sysA_s, sysB_s, "solvent"),
            ]:
                tf_name = f"rbfe_{nameA}_{leg}_{nameB}_{leg}"
                tf = Transformation(
                    stateA=stateA, stateB=stateB,
                    mapping={"ligand": mapping}, protocol=protocol,
                    name=tf_name)

                # Salvage subdir inside the ORIGINAL job directory
                orig_job_dir = prod_dir / str(cluster_id) / tf_name
                salvage_dir = orig_job_dir / "salvage"
                salvage_dir.mkdir(parents=True, exist_ok=True)

                tf.dump(salvage_dir / f"{tf_name}.json")
                shutil.copy(script_src, salvage_dir / "run_quickrun.sh")
                (salvage_dir / "logs").mkdir(exist_ok=True)

                # Provenance sidecar
                prov = {
                    "mapper": "lomap",
                    "element_change": False,
                    "max3d": 1.5,
                    "threed": True,
                    "n_mapped_atoms": n_mapped,
                    "original_mapper": "kartograf",
                    "salvage_reason": "kartograf edge failed",
                    "edge": edge,
                    "leg": leg,
                }
                (salvage_dir / "provenance.json").write_text(
                    json.dumps(prov, indent=2))

                transform_entries.append(f"{cluster_id},{tf_name}")

            n_ok += 1
            print(f"  OK {cluster_id}/{edge}: mapped {n_mapped} atoms")

        except Exception as e:
            print(f"  ERROR {cluster_id}/{edge}: {type(e).__name__}: {e}")
            n_error += 1

    # Append to the queue list so batches accumulate (first batch with
    # batch-start=0 truncates to start fresh; later batches append).
    mode = "w" if args.batch_start == 0 else "a"
    with open(args.outlist, mode) as f:
        if transform_entries:
            f.write("\n".join(transform_entries) + "\n")

    # Append unsalvageable edges similarly
    if skipped_edges:
        import csv
        skip_path = Path(args.outlist).parent / "salvage_skipped_edges.csv"
        write_header = (args.batch_start == 0) or (not skip_path.exists())
        with open(skip_path, "a" if args.batch_start > 0 else "w",
                  newline="") as f:
            w = csv.DictWriter(f, fieldnames=["cluster_id", "edge",
                                              "n_mapped", "frac", "reason"])
            if write_header:
                w.writeheader()
            w.writerows(skipped_edges)
        print(f"Appended {len(skipped_edges)} unsalvageable edges to "
              f"{skip_path}")

    print(f"\n=== Salvage planning batch complete ===")
    print(f"Batch range: {args.batch_start} to "
          f"{args.batch_start + len(batch) - 1}")
    print(f"Edges salvaged this batch (both legs): {n_ok}")
    print(f"No LOMAP mapping found:                 {n_no_mapping}")
    print(f"Low coverage (skipped):                {n_low_coverage}")
    print(f"Errors:                                 {n_error}")
    print(f"Transformation legs written this batch: {len(transform_entries)}")
    next_start = args.batch_start + args.batch_size
    if next_start < total_failed:
        print(f"\nNext batch: --batch-start {next_start}")
    else:
        print(f"\nAll {total_failed} failed edges processed.")
    print(f"Queue list: {args.outlist}")


if __name__ == "__main__":
    main()
