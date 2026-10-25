#!/bin/zsh
# Run this example with Finch V2: 1D and 2D (one Finch session each, under a
# minute in total).
# Usage:  ./run.sh          (from any directory)
# FINCH_V2_DIR is the Finch-V2 checkout (branch mlir-batched-linalg), run with Julia 1.12.
set -e
set -o pipefail          # a Julia failure must not be masked by tee
FINCH_V2_DIR="${FINCH_V2_DIR:-$HOME/Packages/Finch-V2}"
cd "$(dirname "$0")"
: > output/runlog.txt
for d in 1 2; do
    # tee keeps the printed results alongside the other output
    CDI_DIM=$d julia +1.12 --project="$FINCH_V2_DIR" main.jl 2>&1 | tee -a output/runlog.txt
done
