#!/bin/bash
# Run this example with Finch V2 (branch dev).
# Usage:  ./run.sh                          CPU, Julia target, case circle (the reference)
#         CDI_CASE=zalesak ./run.sh         Zalesak's disk; CDI_CASE=sheardrop the vortex drop
#         FINCH_TARGET=mlir_cuda ./run.sh   GPU through the MLIR/IREE target (Linux, IREE set up)
#         FINCH_TARGET=mlir_cpu ./run.sh    the MLIR target on the CPU
# Output lands in output/<case>/ (CDI_OUTDIR overrides), with the printed log in runlog.txt.
# FINCH_V2_DIR is the Finch-V2 checkout (default ~/Packages/Finch-V2), run with Julia 1.12
# (JULIA overrides the command, e.g. JULIA=julia where juliaup's `julia +1.12` is not available).
set -e
set -o pipefail          # a Julia failure must not be masked by tee
FINCH_V2_DIR="${FINCH_V2_DIR:-$HOME/Packages/Finch-V2}"
cd "$(dirname "$0")"
CASE="${CDI_CASE:-circle}"
OUT="${CDI_OUTDIR:-output/$CASE}"
mkdir -p "$OUT"
${JULIA:-julia +1.12} --project="$FINCH_V2_DIR" main.jl 2>&1 | tee "$OUT/runlog.txt"
