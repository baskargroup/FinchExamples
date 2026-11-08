#!/bin/bash
# Submit the MLIR-target runs (Finch's IREE backend) as NOVA jobs, one per
# example, case and target, each ending with the regression check against the
# CPU reference at the solver-tolerance level (rtol 1e-6, atol 1e-9: the MLIR
# target solves with BiCGSTAB + Jacobi).
#
#     bash tools/nova_submit_mlir.sh smoke    Examples 1 and 2 on mlir_cuda and mlir_cpu (4 jobs)
#     bash tools/nova_submit_mlir.sh all      every case on mlir_cuda, plus two mlir_cpu timing runs
#
# mlir_cuda: one A100 (--gres=gpu:a100:1, 8 cores); mlir_cpu: 32 cores (the
# emitted driver runs 32 worker threads). Results go beside the CPU ones:
# output/<case>_<target>/ (single-case examples: output/<target>/; Example 21:
# output_case<X>_<target>/), and the run log in there. Job names
# <ExampleN>_<case>_<target>; names already queued are skipped.
cd "$(dirname "$0")/.."
MODE=${1:-smoke}
submit() {   # folder case target ref
    local folder=$1 case=$2 target=$3 ref=$4
    local short=${folder%%_*}
    local name="${short}${case:+_$case}_$target"
    [ -d "$folder" ] || { echo "skip $name: no folder $folder"; return; }
    squeue -u "$USER" -h -o "%j" | grep -qx "$name" && { echo "skip $name: already queued or running"; return; }
    local outdir
    if [ "$short" = "Example21" ]; then outdir="output_case${case}_$target"
    elif [ -n "$case" ]; then outdir="output/${case}_$target"
    else outdir="output/$target"; fi
    local res
    if [ "$target" = "mlir_cuda" ]; then res="--gres=gpu:a100:1 --cpus-per-task=8 --mem=32G"
    else res="--cpus-per-task=32 --mem=32G"; fi
    sbatch --chdir="$folder" --job-name="$name" --time=04:00:00 $res \
           --export=ALL,FINCH_TARGET=$target,CDI_CASE=$case,CDI_OUTDIR=$outdir,REF=$ref,LOG=$outdir/runlog.txt,RTOL=1e-6,ATOL=1e-9 \
           "$(pwd)/tools/nova_run.sh"
}
E1=Example1_Poisson1D; E2=Example2_Poisson2D; E3=Example3_ConvectionDiffusion1D; E4=Example4_ConvectionDiffusion2D
E5=Example5_ConservativeDiffusiveInterfaceExplicit1D; E6=Example6_ConservativeDiffusiveInterfaceExplicit2D
E7=Example7_ConservativeDiffusiveInterfaceImplicit1D; E8=Example8_ConservativeDiffusiveInterfaceImplicit2D
E21=Example21_NavierStokesConservativeDiffusiveInterfaceExplicit2D
E31=Example31_ConservativeDiffusiveInterfaceScalarTransportExplicit1D; E32=Example32_ConservativeDiffusiveInterfaceScalarTransportExplicit2D
case "$MODE" in
  smoke)
    for t in mlir_cuda mlir_cpu; do submit $E1 "" $t $E1.txt; submit $E2 "" $t $E2.txt; done ;;
  all)
    submit $E1 "" mlir_cuda $E1.txt;  submit $E2 "" mlir_cuda $E2.txt
    submit $E3 "" mlir_cuda $E3.txt;  submit $E4 "" mlir_cuda $E4.txt
    submit $E5 "" mlir_cuda $E5.txt;  submit $E7 "" mlir_cuda $E7.txt
    for c in circle zalesak sheardrop; do submit $E6 $c mlir_cuda ${E6}_$c.txt; submit $E8 $c mlir_cuda ${E8}_$c.txt; done
    submit $E21 A mlir_cuda ${E21}_caseA.txt;  submit $E21 E mlir_cuda ${E21}_caseE.txt
    submit $E31 "" mlir_cuda $E31.txt
    # the 32 references are at 100^2 (CDI_NEL=100): pass it through the environment
    CDI_NEL=100 submit $E32 bubbles mlir_cuda ${E32}_bubbles.txt
    CDI_NEL=100 submit $E32 electrode mlir_cuda ${E32}_electrode.txt
    submit $E6 zalesak mlir_cpu ${E6}_zalesak.txt;  submit $E21 E mlir_cpu ${E21}_caseE.txt ;;
  *) echo "usage: $0 smoke|all"; exit 1 ;;
esac
squeue -u "$USER" -o "%.10i %.28j %.8T %.8M %.14R"
