#!/bin/zsh
# Run this example with the Finch environment.
# Usage:  ./run.sh          (from any directory)
set -e
set -o pipefail          # a Julia failure must not be masked by tee
# The patched Finch that ships with this repository. Set FINCH_DIR in the
# environment to run against a different Finch instead.
FINCH_DIR="${FINCH_DIR:-$(cd "$(dirname "$0")/.." && pwd)/Finch_patched}"
cd "$(dirname "$0")"
# tee keeps the printed results alongside the other output
julia --project="$FINCH_DIR" main.jl 2>&1 | tee output/runlog.txt
