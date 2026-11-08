#!/bin/bash
# Temporal convergence study (see convergence.jl).
# Usage:  ./convergence.sh          (from any directory)
# FINCH_V2_DIR is the Finch-V2 checkout (branch mlir-batched-linalg), run with Julia 1.12
# (JULIA overrides the command, e.g. JULIA=julia where juliaup's `julia +1.12` is not available).
set -e
set -o pipefail          # a Julia failure must not be masked by tee
FINCH_V2_DIR="${FINCH_V2_DIR:-$HOME/Packages/Finch-V2}"
cd "$(dirname "$0")"
# tee keeps the printed results alongside the other output
${JULIA:-julia +1.12} --project="$FINCH_V2_DIR" convergence.jl 2>&1 | tee output/convergence_runlog.txt
