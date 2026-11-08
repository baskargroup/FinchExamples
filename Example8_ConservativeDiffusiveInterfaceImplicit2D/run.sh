#!/bin/bash
# Run this example with Finch V2 (Julia target).
# Usage:  ./run.sh          (from any directory)
# FINCH_V2_DIR is the Finch-V2 checkout (branch mlir-batched-linalg), run with Julia 1.12
# (JULIA overrides the command, e.g. JULIA=julia where juliaup's `julia +1.12` is not available).
set -e
set -o pipefail          # a Julia failure must not be masked by tee
FINCH_V2_DIR="${FINCH_V2_DIR:-$HOME/Packages/Finch-V2}"
cd "$(dirname "$0")"
${JULIA:-julia +1.12} --project="$FINCH_V2_DIR" main.jl 2>&1 | tee output/runlog.txt
