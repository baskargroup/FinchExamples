#!/bin/zsh
# Comparison run: main.jl with the element-wise interface normal instead of
# the nodal one (CDI_NORMAL=elementwise). Results go to
# output/elementwise_normal/.
# Usage:  ./elementwise_normal.sh     (from any directory)
# FINCH_V2_DIR is the Finch-V2 checkout (branch mlir-batched-linalg), run with Julia 1.12.
set -e
set -o pipefail          # a Julia failure must not be masked by tee
FINCH_V2_DIR="${FINCH_V2_DIR:-$HOME/Packages/Finch-V2}"
cd "$(dirname "$0")"
export CDI_NORMAL=elementwise CDI_OUTDIR=output/elementwise_normal
mkdir -p "$CDI_OUTDIR"
# tee keeps the printed results alongside the other output
julia +1.12 --project="$FINCH_V2_DIR" main.jl 2>&1 | tee "$CDI_OUTDIR/runlog.txt"
