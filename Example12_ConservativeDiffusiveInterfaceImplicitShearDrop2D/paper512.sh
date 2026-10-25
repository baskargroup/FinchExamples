#!/bin/zsh
# Drop in shear flow, semi-implicit, on a 512 x 512 mesh with the paper's interface
# parameters (Mirjalili et al. 2020): eps/dx = 0.7, gamma/|u|max = 2.5.
# Results go to output/paper_512/ (several hours).
# Usage:  ./paper512.sh     (from any directory)
# FINCH_V2_DIR is the Finch-V2 checkout (branch mlir-batched-linalg), run with Julia 1.12.
set -e
set -o pipefail          # a Julia failure must not be masked by tee
FINCH_V2_DIR="${FINCH_V2_DIR:-$HOME/Packages/Finch-V2}"
cd "$(dirname "$0")"
export CDI_NEL=512 CDI_EPS_DX=0.7 CDI_GAMMA=2.5 CDI_OUTDIR=output/paper_512
mkdir -p "$CDI_OUTDIR"
# tee keeps the printed results alongside the other output
julia +1.12 --project="$FINCH_V2_DIR" main.jl 2>&1 | tee "$CDI_OUTDIR/runlog.txt"
