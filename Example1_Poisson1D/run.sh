#!/bin/zsh
# Run this example with the Finch environment.
# Usage:  ./run.sh          (from any directory)
set -e
set -o pipefail          # a Julia failure must not be masked by tee
FINCH_DIR="$HOME/Packages/Finch"
cd "$(dirname "$0")"
# tee keeps the printed results alongside the other output
julia --project="$FINCH_DIR" main.jl 2>&1 | tee output/runlog.txt
