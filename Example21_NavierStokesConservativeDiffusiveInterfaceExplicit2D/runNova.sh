#!/bin/bash
# Example 21 on NOVA (SLURM): one case of the 2D Navier-Stokes + CDI run with Finch V2 (Julia target).
# Submit from this folder on NOVA (/work/mech-ai/suresh/projects/FinchExamples/Example21_...):
#
#     sbatch runNova.sh                     case E, output_caseE/
#     sbatch --export=ALL,CDI_CASE=C runNova.sh
#
# Serial Julia run: one task, a few cores for BLAS. Case E took 37 min on the laptop.
# Julia 1.12.7 is in /work/mech-ai/suresh/packages/julia-1.12.7 (official tarball, no module on NOVA); Finch V2 dev is in
# /work/mech-ai/suresh/packages/Finch-V2 (git clone of budvinchathura/Finch-V2, checked out at the commit of the laptop runs).

#SBATCH --time=04:00:00                       # case E about 1 h plus precompilation
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --job-name="ex17_caseE"
#SBATCH --mail-user=sureshmn@iastate.edu
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --output=run_%j.txt                   # job log (stays in this folder)
#SBATCH --error=run_%j.txt

module purge
export PATH=/work/mech-ai/suresh/packages/julia-1.12.7/bin:$PATH
export JULIA_DEPOT_PATH=/work/mech-ai/suresh/packages/julia-depot   # the home quota (10 GB) is too small for the depot
export FINCH_V2_DIR=/work/mech-ai/suresh/packages/Finch-V2
export FINCH_PROJECT=$FINCH_V2_DIR
export JULIA_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}
export CDI_CASE=${CDI_CASE:-E}

cd "$SLURM_SUBMIT_DIR"
pwd; date; julia --version
OUT="output_case${CDI_CASE}"
mkdir -p "$OUT"
# the laptop run.sh is zsh and uses "julia +1.12" (juliaup); same call without them here
julia --project="$FINCH_V2_DIR" main.jl 2>&1 | tee "$OUT/runlog.txt"
echo "exit status ${PIPESTATUS[0]}"
date
