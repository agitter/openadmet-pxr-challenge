#!/bin/bash

set -e

echo "Host: $(hostname)"
echo "Starting OpenMM test..."

python -m openmm.testInstallation
