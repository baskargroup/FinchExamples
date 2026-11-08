#!/bin/bash
# Run this example with Finch V2, development branch (Julia target).
# Usage:  ./run.sh                 case A
#         CDI_CASE=C ./run.sh      case C (A, B, C, D, E: see main.jl)
#         FINCH_TARGET=mlir_cuda ./run.sh   GPU through the MLIR/IREE target (Linux, IREE set up)
# FINCH_V2_DIR is the Finch-V2 checkout on branch dev (default ~/Packages/Finch-V2).
set -e
set -o pipefail          # a Julia failure must not be masked by tee
FINCH_V2_DIR="${FINCH_V2_DIR:-$HOME/Packages/Finch-V2}"
cd "$(dirname "$0")"
CASE="${CDI_CASE:-A}"
OUT="output_case$(echo "$CASE" | tr a-z A-Z)"
[[ -n "$CDI_OUTDIR" ]] && OUT="$CDI_OUTDIR"
mkdir -p "$OUT"
${JULIA:-julia +1.12} --project="$FINCH_V2_DIR" main.jl 2>&1 | tee "$OUT/runlog.txt"
