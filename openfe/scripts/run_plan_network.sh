#!/usr/bin/env bash
#
# openfe/scripts/run_plan_network.sh
#
# HTCondor executable: runs openfe plan-rbfe-network for one cluster.
#
# Arguments:
#   $1 = cluster_id (e.g. "21")

# These two must run BEFORE set -euo pipefail:
#   1. mkdir -p ensures the output dir exists for transfer even on failure
#   2. source activates the conda/mamba env; the activation script
#      references MAMBA_SKIP_ACTIVATE which may be unset, which would
#      cause set -u to abort the script if sourced after set -euo pipefail
mkdir -p network_setup
source /usr/local/bin/_activate_current_env.sh

set -eo pipefail

CLUSTER_ID="$1"

echo "=== Planning network for cluster ${CLUSTER_ID} at $(date) ==="
echo "Host: $(hostname)"
echo "Working directory: $(pwd)"
echo "Scratch contents:"
ls -la

echo "openfe: $(openfe --version)"

LIGANDS="${CLUSTER_ID}/ligands.sdf"
RECEPTOR="${CLUSTER_ID}/receptor.pdb"
SETTINGS="plan_settings.yaml"
OUTDIR="network_setup"

for f in "${LIGANDS}" "${RECEPTOR}" "${SETTINGS}"; do
    if [[ ! -f "${f}" ]]; then
        echo "FATAL: required input not found: ${f}"
        echo "Available files:"
        find . -maxdepth 3 | sort
        exit 1
    fi
done

n_ligands=$(python3 -c "
from rdkit import Chem, RDLogger
RDLogger.DisableLog('rdApp.*')
mols = [m for m in Chem.SDMolSupplier('${LIGANDS}', removeHs=False) if m]
print(len(mols))
")
echo "Ligands in SDF: ${n_ligands}"

openfe plan-rbfe-network \
    -M "${LIGANDS}" \
    -p "${RECEPTOR}" \
    -s "${SETTINGS}" \
    --n-protocol-repeats 1 \
    -n 4 \
    -o "${OUTDIR}" \
    2>&1 | tee "${OUTDIR}/plan.log"

echo "=== Finished cluster ${CLUSTER_ID} at $(date) ==="
echo "Output files in ${OUTDIR}/:"
ls -la "${OUTDIR}/"
