#!/bin/bash
# One example (and case) on NOVA, as a SLURM job, then the regression check.
# Submit from the repository root on NOVA, one job per (example, case):
#
#     sbatch --chdir=Example6_ConservativeDiffusiveInterfaceExplicit2D --job-name=ex6_zalesak \
#            --export=ALL,CDI_CASE=zalesak,REF=Example6_ConservativeDiffusiveInterfaceExplicit2D_zalesak.txt,LOG=output/zalesak/runlog.txt \
#            ../tools/nova_run.sh
#
# tools/nova_submit_all.sh submits every verification run. Variables (all
# through --export): CDI_CASE (empty for single-case examples), REF (file in
# regression/), LOG (the run log relative to the example folder), CDI_NEL etc.
# as the example reads them, FINCH_TARGET (default julia; mlir_cuda needs a
# GPU node: --gres=gpu:a100:1 and the cuda module). The job log nova_<jobid>.txt
# stays in the example folder; it ends with the comparer's verdict.

#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --mail-user=sureshmn@iastate.edu
#SBATCH --mail-type=FAIL
#SBATCH --output=nova_%j.txt
#SBATCH --error=nova_%j.txt

P=/work/mech-ai/suresh/packages
module purge
if [ "${FINCH_TARGET:-julia}" != "julia" ]; then
    module load cmake/3.31.9-onvok7d ninja/1.12.1-py311-l2hkirj cuda/12.8.1-hm53cro
fi
export PATH=$P/julia-1.12.7/bin:$PATH
export JULIA_DEPOT_PATH=$P/julia-depot          # the home quota (10 GB) is too small for the depot
export FINCH_V2_DIR=$P/Finch-V2
export JULIA=julia                              # no juliaup on NOVA: plain julia 1.12.7
export JULIA_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}

pwd; date; hostname; julia --version
echo "example $(basename "$PWD")  case '${CDI_CASE:-}'  target ${FINCH_TARGET:-julia}  nel '${CDI_NEL:-default}'"
start=$(date +%s)
./run.sh > /dev/null 2>&1; status=$?
echo "run.sh exit $status after $(( ($(date +%s) - start) / 60 )) min"
echo "--- summary lines of $LOG"
grep -v -E '^\[|^ *@|^\s*$|Precompil|✓|Finch\]|^ [├└─ ]|^─|Finch has completed|dependenc|precompiled' "$LOG" | tail -20
if [ -n "$REF" ]; then
    echo "--- regression check against regression/$REF (rtol ${RTOL:-1e-8}, atol ${ATOL:-1e-12})"
    python3 ../tools/compare_runlog.py --rtol "${RTOL:-1e-8}" --atol "${ATOL:-1e-12}" "../regression/$REF" "$LOG"
    echo "compare exit $?"
fi
date
