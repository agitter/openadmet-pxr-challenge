#!/usr/bin/env bash
#
# openfe/scripts/run_plan_network.sh
#
# HTCondor executable: runs openfe plan-rbfe-network for one cluster.
#
# Arguments:
#   $1 = cluster_id (e.g. "21")
#
# HTCondor transfers plan_inputs/<cluster_id>/ as a flat directory
# named <cluster_id>/ in the scratch dir (parent path stripped).
# plan_settings.yaml lands at the top level alongside it.
#
# Writes output to ./network_setup/ in scratch; transfer_output_remaps
# in the submit file redirects this to network_setup/<cluster_id>/ on
# the access point.
#
# Runs inside: osdf:///chtc/staging/a/agitter/containers/openfe_1.11.1.sif

set -euo pipefail

CLUSTER_ID="$1"

echo "=== Planning network for cluster ${CLUSTER_ID} at $(date) ==="
echo "Host: $(hostname)"
echo "Working directory: $(pwd)"
echo "Scratch contents:"
ls -la

source /usr/local/bin/_activate_current_env.sh

echo "openfe: $(openfe --version)"

# Input paths - HTCondor strips the parent dir on transfer so
# plan_inputs/<cluster_id>/ arrives as <cluster_id>/ here.
LIGANDS="${CLUSTER_ID}/ligands.sdf"
RECEPTOR="${CLUSTER_ID}/receptor.pdb"
SETTINGS="plan_settings.yaml"

# Output dir in scratch - remapped to network_setup/<cluster_id>/
# on the access point via transfer_output_remaps.
OUTDIR="network_setup"

# Pre-create output dir so transfer never fails with "no such file"
# even if openfe itself errors out.
mkdir -p "${OUTDIR}"

# Verify inputs exist
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
    --verbose \
    2>&1 | tee "${OUTDIR}/plan.log"

echo "=== Finished cluster ${CLUSTER_ID} at $(date) ==="
echo "Output files in ${OUTDIR}/:"
ls -la "${OUTDIR}/"
