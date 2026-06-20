#!/usr/bin/env bash
#
# openfe/scripts/run_quickrun_timing.sh
#
# Timing test: runs openfe quickrun on a single transformation JSON,
# measuring wall-clock time per leg. Used to calibrate GPUJobLength
# and sampling settings before launching the full 1066-job campaign.
#
# Arguments:
#   $1 = transformation JSON filename (e.g. "rbfe_OCNT-2310728_complex_OADMET-0006495_complex.json")
#        Note: openfe orders ligands in the filename based on MST edge direction,
#        so the anchor compound typically appears first.
#
# Files expected in scratch (from transfer_input_files):
#   $(transformation)      the transformation JSON
#   network_setup.json     the alchemical network JSON (needed by quickrun)
#
# Output:
#   quickrun_output/       openfe quickrun results directory
#   timing.txt             wall-clock timing summary

mkdir -p quickrun_output
# Create timing.txt early so it always exists for transfer, even on eviction.
# It will be overwritten with real content on successful completion.
touch timing.txt
source /usr/local/bin/_activate_current_env.sh

set -eo pipefail

TRANSFORMATION="$1"
echo "=== openfe quickrun timing test ==="
echo "Host: $(hostname)"
echo "Transformation: ${TRANSFORMATION}"
echo "Start: $(date)"
START=$(date +%s)

# Report GPU info
nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader || echo "nvidia-smi not available"

echo "openfe: $(openfe --version)"
echo "Python: $(python --version)"

# Report OpenMM platform
python -c "
import openmm
print('OpenMM version:', openmm.__version__)
from openmm import Platform
for i in range(Platform.getNumPlatforms()):
    p = Platform.getPlatform(i)
    print(f'  Platform {i}: {p.getName()}')
"

openfe quickrun \
    "${TRANSFORMATION}" \
    -d quickrun_output \
    -o quickrun_output/result.json \
    --resume

END=$(date +%s)
ELAPSED=$((END - START))

echo ""
echo "=== Timing Summary ==="
echo "Transformation: ${TRANSFORMATION}"
echo "Wall-clock time: ${ELAPSED} seconds ($(( ELAPSED / 60 )) min $(( ELAPSED % 60 )) sec)"
echo "Host: $(hostname)"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader 2>/dev/null || true
echo "End: $(date)"

# Write timing file for easy inspection
cat > timing.txt << TIMING
transformation: ${TRANSFORMATION}
wall_clock_seconds: ${ELAPSED}
wall_clock_human: $(( ELAPSED / 60 ))m $(( ELAPSED % 60 ))s
host: $(hostname)
gpu: $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null || echo unknown)
TIMING

echo "Wrote timing.txt"
ls -la quickrun_output/
