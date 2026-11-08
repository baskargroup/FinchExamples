#!/bin/bash
# Run this example with Finch V2 (branch dev).
#         FINCH_TARGET=mlir_cuda ./run.sh   GPU through the MLIR/IREE target (Linux, IREE set up)
#         FINCH_TARGET=mlir_cpu ./run.sh    the MLIR target on the CPU
# Usage:  ./run.sh          (from any directory)
# FINCH_V2_DIR is the Finch-V2 checkout (default ~/Packages/Finch-V2), run with Julia 1.12
# (JULIA overrides the command, e.g. JULIA=julia where juliaup's `julia +1.12` is not available).
set -e
set -o pipefail          # a Julia failure must not be masked by tee
FINCH_V2_DIR="${FINCH_V2_DIR:-$HOME/Packages/Finch-V2}"
cd "$(dirname "$0")"
OUT="${CDI_OUTDIR:-output}"            # CDI_OUTDIR: another results folder (e.g. output/mlir_cuda)
mkdir -p "$OUT"
${JULIA:-julia +1.12} --project="$FINCH_V2_DIR" main.jl 2>&1 | tee "$OUT/runlog.txt"
