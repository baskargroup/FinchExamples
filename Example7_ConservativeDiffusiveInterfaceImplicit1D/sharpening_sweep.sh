#!/bin/bash
# Sharpening-term study (see sharpening_sweep.jl): both meshes, one after the other.
# Usage:  ./sharpening_sweep.sh          (from any directory)
# FINCH_V2_DIR is the Finch-V2 checkout (branch mlir-batched-linalg), run with Julia 1.12
# (JULIA overrides the command, e.g. JULIA=julia where juliaup's `julia +1.12` is not available).
set -e
set -o pipefail          # a Julia failure must not be masked by tee
FINCH_V2_DIR="${FINCH_V2_DIR:-$HOME/Packages/Finch-V2}"
cd "$(dirname "$0")"
mkdir -p output/sharpening_sweep
for nel in 100 200; do
    ${JULIA:-julia +1.12} --project="$FINCH_V2_DIR" sharpening_sweep.jl $nel 2>&1 | tee -a output/sharpening_sweep/runlog.txt
done
