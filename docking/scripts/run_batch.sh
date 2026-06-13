#!/usr/bin/env bash
#
# docking/scripts/run_batch.sh  (v3 - thin orchestrator, no embedded
# Python; calls embed_ligand.py and dock_all_receptors.py, both
# transferred alongside this script as separate input files)
#
# HTCondor executable for one docking job. Each job receives:
#   $1 = work unit file (e.g. work_units/lig_0259.csv)
#         single-row CSV: cluster_id,ligand_name,ligand_smiles
#
# Other files expected in the job's working directory (all listed in
# transfer_input_files in submit_docking.sub):
#   - receptors.zip          (all ~62 receptor PDBQTs + boxes.csv)
#   - embed_ligand.py         (SMILES -> 3D SDF)
#   - dock_all_receptors.py   (loops over receptors, runs gnina)
#
# Expected runtime: ~62 receptors x 1-5 min/docking = ~62-310 min/job.
#
# Output: results/<cluster_id>/<cluster_id>__<pdb_id>_docked.sdf.gz
#         results/<cluster_id>/<cluster_id>__<pdb_id>.log(.stdout)
#         results/<cluster_id>_summary.csv
#
# Runs INSIDE gnina/gnina:v1.3.1 (has gnina + obabel + python3).

set -euo pipefail

WORK_UNIT="$1"

# Parse the single data row (skip header)
IFS=, read -r cluster_id ligand_name ligand_smiles < <(tail -n +2 "$WORK_UNIT")

echo "=== Starting ligand ${ligand_name} (cluster ${cluster_id}) at $(date) ==="

mkdir -p "results/${cluster_id}"

# ---------------------------------------------------------------
# Ensure python deps available (rdkit + meeko + pandas for ligand prep).
# The gnina/gnina:v1.3.1 image ships Python 3.8.10.
# ---------------------------------------------------------------
python3 -c "import rdkit, meeko, pandas" 2>/dev/null || {
    echo "Installing rdkit + meeko + pandas..."
    pip install --user --quiet "rdkit==2023.9.6" "meeko==0.5.0" "pandas==2.0.3"
}
python3 -c "import rdkit, meeko, pandas" || {
    echo "FATAL: required python deps unavailable after install attempt"
    exit 1
}

# --user installs place console scripts (mk_prepare_ligand.py etc.) in
# ~/.local/bin, which is not on PATH by default.
export PATH="${HOME}/.local/bin:${PATH}"

# Log the versions actually in use for this job
python3 -c "
import rdkit, meeko, pandas
print(f'rdkit={rdkit.__version__}')
print(f'meeko={meeko.__version__}')
print(f'pandas={pandas.__version__}')
" > "results/${cluster_id}/dep_versions.txt"
cat "results/${cluster_id}/dep_versions.txt"

# ---------------------------------------------------------------
# Unzip the shared receptor archive (one file transferred per job)
# ---------------------------------------------------------------
mkdir -p receptors
python3 -c "import zipfile; zipfile.ZipFile('receptors.zip').extractall('receptors')"
n_receptors=$(find receptors -name "*_protein.pdbqt" | wc -l)
echo "Unzipped ${n_receptors} receptor PDBQTs"

# ---------------------------------------------------------------
# 1. SMILES -> 3D SDF -> PDBQT (once per ligand)
# ---------------------------------------------------------------
mkdir -p ligand_3d
lig_sdf="ligand_3d/${cluster_id}.sdf"
lig_pdbqt="ligand_3d/${cluster_id}.pdbqt"

python3 embed_ligand.py "$ligand_smiles" "$lig_sdf"

mk_prepare_ligand.py -i "$lig_sdf" -o "$lig_pdbqt" \
    > "results/${cluster_id}/ligand_prep.log" 2>&1 || {
    echo "FATAL: ligand prep failed for cluster ${cluster_id}, "
         "see results/${cluster_id}/ligand_prep.log"
    exit 1
}

# ---------------------------------------------------------------
# 2. Dock against every receptor (loop + gnina calls handled in Python)
# ---------------------------------------------------------------
python3 dock_all_receptors.py \
    --ligand-pdbqt "$lig_pdbqt" \
    --boxes receptors/boxes.csv \
    --receptor-dir receptors \
    --cluster-id "$cluster_id" \
    --ligand-name "$ligand_name" \
    --outdir "results/${cluster_id}" \
    --summary "results/${cluster_id}_summary.csv" \
    --exhaustiveness 8 \
    --num-modes 5

# clean up scratch to keep transferred output small
rm -rf receptors ligand_3d

echo "=== Finished ligand ${ligand_name} (cluster ${cluster_id}) at $(date) ==="
