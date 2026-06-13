#!/usr/bin/env bash
#
# docking/scripts/run_batch.sh
#
# HTCondor executable for one docking batch job. Each job receives:
#   $1 = work unit file (e.g. work_units/batch_0042.csv)
#
# The work unit CSV has columns: cluster_id,ligand_name,ligand_smiles,pdb_id
# i.e. one row per (cluster representative, receptor) docking task.
#
# For each row, this script:
#   1. Converts the ligand SMILES to a 3D SDF with RDKit (ETKDG embedding
#      + MMFF optimization), then to PDBQT with Meeko.
#   2. Looks up the receptor PDBQT and docking box (center/size) from
#      receptors/boxes.csv (staged alongside this job).
#   3. Runs gnina with the CNN scoring enabled, writing the top poses
#      (SDF) and a score summary.
#
# All (ligand, receptor) pairs in the work unit run sequentially within
# this single job - this is what gives each job a ~20-100 minute runtime
# instead of many tiny sub-minute jobs.
#
# Output: results/<batch_name>/<cluster_id>__<pdb_id>.sdf.gz
#         results/<batch_name>/<cluster_id>__<pdb_id>.log
#         results/<batch_name>_summary.csv
#
# Runs INSIDE the gnina/gnina:v1.3.1 container (has gnina + obabel +
# python3). Meeko/RDKit are installed at job start via pip if not
# already present in the image.

set -euo pipefail

WORK_UNIT="$1"
BATCH_NAME=$(basename "${WORK_UNIT}" .csv)

echo "=== Starting batch ${BATCH_NAME} at $(date) ==="

mkdir -p "results/${BATCH_NAME}"

# ---------------------------------------------------------------
# Ensure python deps available (gemmi + rdkit + meeko for ligand prep).
# The gnina/gnina image ships its own python; pip install --user
# to avoid permission issues.
# ---------------------------------------------------------------
python3 -c "import gemmi, rdkit, meeko" 2>/dev/null || {
    echo "Installing gemmi + rdkit + meeko ..."
    pip install --user --quiet gemmi rdkit meeko
}
python3 -c "import gemmi, rdkit, meeko" || {
    echo "FATAL: required python deps unavailable after install attempt"
    exit 1
}

# ---------------------------------------------------------------
# Per-pair docking function
# ---------------------------------------------------------------
dock_one() {
    local cluster_id="$1"
    local ligand_name="$2"
    local ligand_smiles="$3"
    local pdb_id="$4"

    local tag="${cluster_id}__${pdb_id}"
    local lig_sdf="ligands_3d/${tag}.sdf"
    local lig_pdbqt="ligands_3d/${tag}.pdbqt"
    local out_sdf="results/${BATCH_NAME}/${tag}_docked.sdf"
    local out_log="results/${BATCH_NAME}/${tag}.log"

    mkdir -p ligands_3d

    # 1. SMILES -> 3D SDF (RDKit ETKDG + MMFF)
    python3 - "$ligand_smiles" "$lig_sdf" <<'PYEOF'
import sys
from rdkit import Chem
from rdkit.Chem import AllChem

smi, out_path = sys.argv[1], sys.argv[2]
mol = Chem.MolFromSmiles(smi)
if mol is None:
    sys.exit(f"Failed to parse SMILES: {smi}")
mol = Chem.AddHs(mol)
params = AllChem.ETKDGv3()
params.randomSeed = 42
cid = AllChem.EmbedMolecule(mol, params)
if cid < 0:
    sys.exit(f"Embedding failed for SMILES: {smi}")
try:
    AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
except Exception:
    AllChem.UFFOptimizeMolecule(mol, maxIters=500)

w = Chem.SDWriter(out_path)
w.write(mol)
w.close()
PYEOF

    # 2. SDF -> PDBQT (Meeko)
    mk_prepare_ligand.py -i "$lig_sdf" -o "$lig_pdbqt" \
        > "${out_log}.ligprep" 2>&1 || {
        echo "  [${tag}] ligand prep FAILED, see ${out_log}.ligprep"
        return 1
    }

    # 3. Look up receptor + box
    local receptor="receptors/${pdb_id}_protein.pdbqt"
    if [[ ! -f "$receptor" ]]; then
        echo "  [${tag}] receptor not found: ${receptor}"
        return 1
    fi

    read -r cx cy cz sx sy sz <<< "$(python3 - "$pdb_id" <<'PYEOF'
import sys, csv
pdb_id = sys.argv[1]
with open("receptors/boxes.csv") as f:
    for row in csv.DictReader(f):
        if row["pdb_id"] == pdb_id:
            print(row["center_x"], row["center_y"], row["center_z"],
                  row["size_x"], row["size_y"], row["size_z"])
            break
PYEOF
)"

    # 4. Dock with GNINA (CNN scoring enabled, default exhaustiveness)
    gnina \
        -r "$receptor" \
        -l "$lig_pdbqt" \
        --center_x "$cx" --center_y "$cy" --center_z "$cz" \
        --size_x "$sx" --size_y "$sy" --size_z "$sz" \
        --exhaustiveness 8 \
        --num_modes 5 \
        --cnn_scoring rescore \
        -o "$out_sdf" \
        --log "$out_log" \
        > "${out_log}.stdout" 2>&1 || {
        echo "  [${tag}] gnina FAILED, see ${out_log}.stdout"
        return 1
    }

    gzip -f "$out_sdf" 2>/dev/null || true

    echo "  [${tag}] done"
    return 0
}

# ---------------------------------------------------------------
# Iterate over work unit rows (skip header)
# ---------------------------------------------------------------
SUMMARY="results/${BATCH_NAME}_summary.csv"
echo "cluster_id,ligand_name,pdb_id,status" > "$SUMMARY"

tail -n +2 "$WORK_UNIT" | while IFS=, read -r cluster_id ligand_name ligand_smiles pdb_id; do
    echo "--- Docking cluster ${cluster_id} (${ligand_name}) into ${pdb_id} ---"
    if dock_one "$cluster_id" "$ligand_name" "$ligand_smiles" "$pdb_id"; then
        echo "${cluster_id},${ligand_name},${pdb_id},success" >> "$SUMMARY"
    else
        echo "${cluster_id},${ligand_name},${pdb_id},FAILED" >> "$SUMMARY"
    fi
done

echo "=== Finished batch ${BATCH_NAME} at $(date) ==="
echo "Summary written to ${SUMMARY}"
