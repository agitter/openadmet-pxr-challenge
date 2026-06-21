#!/usr/bin/env bash
#
# openfe/scripts/run_quickrun.sh
#
# HTCondor executable for production openfe quickrun jobs.
# No checkpointing - jobs are short enough (~1.6h mean, 7h max) to
# complete in a single GPU slot within the 12h "short" time limit.
# Simpler than checkpoint/resume and avoids the -o file conflict bug.
#
# Arguments:
#   $1 = transformation JSON filename
#
# Files transferred from initialdir:
#   $1                  the transformation JSON
#   network_setup.json  the alchemical network JSON
#
# Output transferred back to initialdir:
#   result.json         the dG estimate (~3 MB on success, 0 bytes on failure)

# Create result.json early so HTCondor transfer never fails on a missing file
touch result.json

# Activate conda/mamba environment before set -eo pipefail
source /usr/local/bin/_activate_current_env.sh

set -eo pipefail

TRANSFORMATION="$1"
echo "=== openfe quickrun production ==="
echo "Host: $(hostname)"
echo "Transformation: ${TRANSFORMATION}"
echo "Start: $(date)"
START=$(date +%s)

nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader 2>/dev/null || echo "nvidia-smi not available"
echo "openfe: $(openfe --version)"

# CUDA fast-fail: exit immediately if CUDA not available.
# Non-zero exit triggers HTCondor retry on a different server.
python -c "
import openmm
from openmm import Platform
platforms = [Platform.getPlatform(i).getName() for i in range(Platform.getNumPlatforms())]
print('Available platforms:', platforms)
if 'CUDA' not in platforms:
    raise RuntimeError('CUDA not available: ' + str(platforms))
print('CUDA confirmed available.')
"

# Run without --resume: fresh start every time, no checkpoint files.
# Work directory is a local scratch dir that stays on the execute node.
mkdir -p quickrun_scratch
openfe quickrun \
    "${TRANSFORMATION}" \
    -d quickrun_scratch \
    -o result.json

END=$(date +%s)
ELAPSED=$((END - START))
echo ""
echo "=== Completed ==="
echo "Wall-clock: ${ELAPSED}s ($(( ELAPSED / 3600 ))h $(( (ELAPSED % 3600) / 60 ))m)"
echo "result.json size: $(wc -c < result.json) bytes"
echo "End: $(date)"
