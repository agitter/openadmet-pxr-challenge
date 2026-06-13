#!/bin/bash

set -eo pipefail

source /usr/local/bin/_activate_current_env.sh

echo "Host: $(hostname)"

echo "Python: $(command -v python)"
python --version

echo "OpenFE: $(command -v openfe)"
openfe --version

python -m openmm.testInstallation
