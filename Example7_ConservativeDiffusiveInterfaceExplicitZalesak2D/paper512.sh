#!/bin/zsh
# Zalesak's slotted disk, explicit, on a 512 x 512 mesh with the interface
# parameters of the drop-in-shear 512 runs (Mirjalili et al. 2020): eps/dx = 0.7,
# gamma/|u|max = 2.5. Results go to output/paper_512/ (about 4-5 hours).
# Usage:  ./paper512.sh     (from any directory)
set -e
set -o pipefail          # a Julia failure must not be masked by tee
# The patched Finch that ships with this repository. Set FINCH_DIR in the
# environment to run against a different Finch instead.
FINCH_DIR="${FINCH_DIR:-$(cd "$(dirname "$0")/.." && pwd)/Finch_patched}"
cd "$(dirname "$0")"
export CDI_NEL=512 CDI_EPS_DX=0.7 CDI_GAMMA=2.5 CDI_OUTDIR=output/paper_512
mkdir -p "$CDI_OUTDIR"
# tee keeps the printed results alongside the other output
julia --project="$FINCH_DIR" main.jl 2>&1 | tee "$CDI_OUTDIR/runlog.txt"
