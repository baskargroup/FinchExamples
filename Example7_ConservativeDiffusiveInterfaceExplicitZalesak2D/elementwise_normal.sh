#!/bin/zsh
# Comparison run: main.jl with the element-wise interface normal instead of
# the nodal one (see elementwise_normal.jl). Results go to
# output/elementwise_normal/ (about 10 minutes on a laptop).
# Usage:  ./elementwise_normal.sh     (from any directory)
set -e
set -o pipefail          # a Julia failure must not be masked by tee
# The patched Finch that ships with this repository. Set FINCH_DIR in the
# environment to run against a different Finch instead.
FINCH_DIR="${FINCH_DIR:-$(cd "$(dirname "$0")/.." && pwd)/Finch_patched}"
cd "$(dirname "$0")"
mkdir -p output/elementwise_normal
# tee keeps the printed results alongside the other output
julia --project="$FINCH_DIR" elementwise_normal.jl 2>&1 | tee output/elementwise_normal/runlog.txt
