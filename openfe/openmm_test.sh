#!/bin/bash

set -euo pipefail

nvidia-smi || true

source /usr/local/bin/_activate_current_env.sh

echo "Python:"
command -v python
python --version

echo "OpenFE:"
command -v openfe
openfe --version

python -m openmm.testInstallation
