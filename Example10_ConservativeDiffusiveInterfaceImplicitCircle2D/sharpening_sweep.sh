#!/bin/zsh
# Sharpening-term study (see sharpening_sweep.jl): both meshes, one after the other.
# Usage:  ./sharpening_sweep.sh          (from any directory)
set -e
set -o pipefail          # a Julia failure must not be masked by tee
# The patched Finch that ships with this repository. Set FINCH_DIR in the
# environment to run against a different Finch instead.
FINCH_DIR="${FINCH_DIR:-$(cd "$(dirname "$0")/.." && pwd)/Finch_patched}"
cd "$(dirname "$0")"
mkdir -p output/sharpening_sweep
for nel in 100 200; do
    julia --project="$FINCH_DIR" sharpening_sweep.jl $nel 2>&1 | tee -a output/sharpening_sweep/runlog.txt
done
