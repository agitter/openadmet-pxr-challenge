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

Usage (inside container, on AP):
    python openfe/scripts/10_salvage_plan.py \
        --edge-report openfe/final_edge_report.csv \
        --production-dir openfe/production \
        --rbfe-inputs openfe/rbfe_inputs \
        --script openfe/scripts/run_quickrun.sh

All edges with both_done=False are salvaged (every incomplete edge).
After the salvage run completes, re-run connectivity analysis on the
combined Kartograf + LOMAP results to see total compound coverage.
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
    ap.add_argument("--outlist", default="openfe/production/salvage_transform_list.txt",
                    help="Output queue list for the salvage production run")
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
    salvage = salvage[~salvage["both_done"]].copy()
    print(f"Salvaging {len(salvage)} failed edges (all incomplete edges)")

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
    n_ok = n_no_mapping = n_error = 0

    for _, row in salvage.iterrows():
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

    with open(args.outlist, "w") as f:
        f.write("\n".join(transform_entries) + "\n")

    print(f"\n=== Salvage planning complete ===")
    print(f"Edges salvaged (both legs):  {n_ok}")
    print(f"No LOMAP mapping found:      {n_no_mapping}")
    print(f"Errors:                      {n_error}")
    print(f"Transformation legs written: {len(transform_entries)}")
    print(f"Wrote queue list: {args.outlist}")


if __name__ == "__main__":
    main()
