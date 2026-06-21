#!/usr/bin/env bash
#
# openfe/scripts/run_quickrun.sh
#
# HTCondor executable for production openfe quickrun jobs.
# Runs a single RBFE transformation leg (complex or solvent).
#
# Arguments:
#   $1 = transformation JSON filename (e.g. "rbfe_OCNT-2310728_complex_OADMET-0006495_complex.json")
#
# Files expected in scratch (transferred from initialdir):
#   $1                     the transformation JSON
#   network_setup.json     the alchemical network JSON
#   quickrun_output/       empty on first run, contains checkpoints on resume
#
# Outputs (transferred back to initialdir):
#   quickrun_output/       checkpoint files during run, cleaned after success
#   result.json            the ΔG estimate (copied out of quickrun_output/ on success)
#   COMPLETED              marker file indicating successful completion
#   FAILED                 marker file with error info if job fails permanently

# Pre-failure-safe setup: create output files/dirs before anything can fail,
# so HTCondor transfer never errors on missing paths.
mkdir -p quickrun_output
touch result.json

# Activate conda/mamba environment (must be before set -eo pipefail
# because the activation script references unset variables)
source /usr/local/bin/_activate_current_env.sh

set -eo pipefail

TRANSFORMATION="$1"
echo "=== openfe quickrun production ==="
echo "Host: $(hostname)"
echo "Transformation: ${TRANSFORMATION}"
echo "Start: $(date)"
START=$(date +%s)

# GPU info (nvidia-smi may not be available in all containers)
nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader 2>/dev/null || echo "nvidia-smi not available"

echo "openfe: $(openfe --version)"

# CUDA fast-fail check: exit immediately if CUDA not available rather
# than running for hours on CPU and producing NaN errors.
python -c "
import openmm
from openmm import Platform
platforms = [Platform.getPlatform(i).getName() for i in range(Platform.getNumPlatforms())]
print('Available platforms:', platforms)
if 'CUDA' not in platforms:
    raise RuntimeError('CUDA platform not available - job must run on a GPU. '
                       'Available platforms: ' + str(platforms))
print('CUDA platform confirmed available.')
"

# Check if this is a resume (quickrun_cache exists from prior run)
if [ -d "quickrun_output/quickrun_cache" ]; then
    echo "Found existing quickrun_cache - this is a RESUME run"
else
    echo "No quickrun_cache found - this is a FRESH run"
fi

# Run the transformation with --resume (safe even on fresh runs)
openfe quickrun \
    "${TRANSFORMATION}" \
    -d quickrun_output \
    -o quickrun_output/result.json \
    --resume

END=$(date +%s)
ELAPSED=$((END - START))

echo ""
echo "=== Completed successfully ==="
echo "Wall-clock: ${ELAPSED} seconds ($(( ELAPSED / 3600 ))h $(( (ELAPSED % 3600) / 60 ))m)"
echo "End: $(date)"

# Copy result.json to top level (alongside quickrun_output/ for transfer)
cp quickrun_output/result.json result.json

# Clean up large simulation files to save AP storage.
# Keep only result.json and the quickrun_cache (small, needed if
# HTCondor somehow restarts this job after completion).
# The large files are in shared_*SimulationUnit*/ directories.
echo "Cleaning up large checkpoint/trajectory files..."
du -sh quickrun_output/ 2>/dev/null || true
find quickrun_output -type d -name "shared_*SimulationUnit*" -exec rm -rf {} + 2>/dev/null || true
find quickrun_output -name "*.nc" -delete 2>/dev/null || true
du -sh quickrun_output/ 2>/dev/null || true

# Write completion marker
echo "transformation: ${TRANSFORMATION}" > COMPLETED
echo "wall_clock_seconds: ${ELAPSED}" >> COMPLETED
echo "host: $(hostname)" >> COMPLETED
echo "gpu: $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null || echo unknown)" >> COMPLETED
echo "end_time: $(date)" >> COMPLETED
