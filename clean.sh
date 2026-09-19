#!/bin/zsh
# Delete simulation outputs: every file that git ignores (CSV, logs, PNG, MP4, ...).
# Files tracked by git or not covered by a .gitignore are never touched.
#
# Usage:
#   ./clean.sh                                   # clean all example folders
#   ./clean.sh Example3_ConvectionDiffusion1D    # clean one folder
#   ./clean.sh -n                                # dry run: list what would be deleted
#   ./clean.sh -n Example3_ConvectionDiffusion1D
set -e
cd "$(dirname "$0")"
git clean -fdX "$@"
