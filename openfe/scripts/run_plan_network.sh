#!/usr/bin/env bash
#
# openfe/scripts/run_plan_network.sh
#
# HTCondor executable: runs openfe plan-rbfe-network for one cluster.
#
# Each job receives:
#   $1 = cluster_id (e.g. "21")
#
# Files expected in the job's working directory
# (listed in transfer_input_files in submit_plan_network.sub):
#   plan_inputs/<cluster_id>/ligands.sdf
#   plan_inputs/<cluster_id>/receptor.pdb
#   plan_settings.yaml
#
# Output:
#   network_setup/<cluster_id>/          openfe network JSON files
#   network_setup/<cluster_id>.log       planning log
#
# Runs INSIDE the openfe Singularity image (openfe_1.11.1.sif).

set -euo pipefail

CLUSTER_ID="$1"

echo "=== Planning network for cluster ${CLUSTER_ID} at $(date) ==="
echo "Host: $(hostname)"

source /usr/local/bin/_activate_current_env.sh

echo "openfe version: $(openfe --version)"

LIGANDS="plan_inputs/${CLUSTER_ID}/ligands.sdf"
RECEPTOR="plan_inputs/${CLUSTER_ID}/receptor.pdb"
SETTINGS="plan_settings.yaml"
OUTDIR="network_setup/${CLUSTER_ID}"

# Verify inputs
for f in "$LIGANDS" "$RECEPTOR" "$SETTINGS"; do
    if [[ ! -f "$f" ]]; then
        echo "FATAL: required input not found: $f"
        exit 1
    fi
done

n_ligands=$(python3 -c "
from rdkit import Chem
mols = [m for m in Chem.SDMolSupplier('$LIGANDS', removeHs=False) if m]
print(len(mols))
")
echo "Ligands in SDF: ${n_ligands}"

mkdir -p "$OUTDIR"

openfe plan-rbfe-network \
    -M "$LIGANDS" \
    -p "$RECEPTOR" \
    -s "$SETTINGS" \
    --n-protocol-repeats 1 \
    -n 4 \
    -o "$OUTDIR" \
    --verbose \
    2>&1 | tee "network_setup/${CLUSTER_ID}.log"

echo "=== Finished cluster ${CLUSTER_ID} at $(date) ==="
echo "Output files:"
ls -la "$OUTDIR/"
