#!/bin/zsh
# Run this example with the Finch environment.
# Usage:  ./run.sh          (from any directory)
set -e
FINCH_DIR="$HOME/Packages/Finch"
cd "$(dirname "$0")"
julia --project="$FINCH_DIR" main.jl
