#!/bin/bash
# Submit every verification run as its own NOVA job (tools/nova_run.sh).
# Run from the repository root on NOVA after `git pull`:   bash tools/nova_submit_all.sh
# Each line: example folder | case | extra env | reference file | run log | rtol | atol
cd "$(dirname "$0")/.."
submit() {   # folder case extra ref log rtol atol
    local folder=$1 case=$2 extra=$3 ref=$4 log=$5 rtol=${6:-1e-8} atol=${7:-1e-12}
    local name="${folder%%_*}${case:+_$case}"
    [ -d "$folder" ] || { echo "skip $name: no folder $folder (not on this branch yet)"; return; }
    sbatch --chdir="$folder" --job-name="$name" \
           --export=ALL,CDI_CASE=$case,REF=$ref,LOG=$log,RTOL=$rtol,ATOL=$atol${extra:+,$extra} \
           ../tools/nova_run.sh
}
submit Example1_Poisson1D                                            ""        "" Example1_Poisson1D.txt                                               output/runlog.txt
submit Example2_Poisson2D                                            ""        "" Example2_Poisson2D.txt                                               output/runlog.txt
submit Example3_ConvectionDiffusion1D                                ""        "" Example3_ConvectionDiffusion1D.txt                                   output/runlog.txt
submit Example4_ConvectionDiffusion2D                                ""        "" Example4_ConvectionDiffusion2D.txt                                   output/runlog.txt
submit Example5_ConservativeDiffusiveInterfaceExplicit1D             ""        "" Example5_ConservativeDiffusiveInterfaceExplicit1D.txt                output/runlog.txt
submit Example6_ConservativeDiffusiveInterfaceExplicit2D             circle    "" Example6_ConservativeDiffusiveInterfaceExplicit2D_circle.txt         output/circle/runlog.txt
submit Example6_ConservativeDiffusiveInterfaceExplicit2D             zalesak   "" Example6_ConservativeDiffusiveInterfaceExplicit2D_zalesak.txt        output/zalesak/runlog.txt
submit Example6_ConservativeDiffusiveInterfaceExplicit2D             sheardrop "" Example6_ConservativeDiffusiveInterfaceExplicit2D_sheardrop.txt      output/sheardrop/runlog.txt
submit Example7_ConservativeDiffusiveInterfaceImplicit1D             ""        "" Example7_ConservativeDiffusiveInterfaceImplicit1D.txt                output/runlog.txt
submit Example8_ConservativeDiffusiveInterfaceImplicit2D             circle    "" Example8_ConservativeDiffusiveInterfaceImplicit2D_circle.txt         output/circle/runlog.txt
submit Example8_ConservativeDiffusiveInterfaceImplicit2D             zalesak   "" Example8_ConservativeDiffusiveInterfaceImplicit2D_zalesak.txt        output/zalesak/runlog.txt
submit Example8_ConservativeDiffusiveInterfaceImplicit2D             sheardrop "" Example8_ConservativeDiffusiveInterfaceImplicit2D_sheardrop.txt      output/sheardrop/runlog.txt
submit Example21_NavierStokesConservativeDiffusiveInterfaceExplicit2D A         "" Example21_NavierStokesConservativeDiffusiveInterfaceExplicit2D_caseA.txt output_caseA/runlog.txt
submit Example21_NavierStokesConservativeDiffusiveInterfaceExplicit2D E         "" Example21_NavierStokesConservativeDiffusiveInterfaceExplicit2D_caseE.txt output_caseE/runlog.txt
submit Example31_ConservativeDiffusiveInterfaceScalarTransportExplicit1D ""     "" Example31_ConservativeDiffusiveInterfaceScalarTransportExplicit1D.txt output/runlog.txt 1e-6 1e-9
submit Example32_ConservativeDiffusiveInterfaceScalarTransportExplicit2D bubbles   CDI_NEL=100 Example32_ConservativeDiffusiveInterfaceScalarTransportExplicit2D_bubbles.txt   output/bubbles/runlog.txt   1e-6 1e-9
submit Example32_ConservativeDiffusiveInterfaceScalarTransportExplicit2D electrode CDI_NEL=100 Example32_ConservativeDiffusiveInterfaceScalarTransportExplicit2D_electrode.txt output/electrode/runlog.txt 1e-6 1e-9
squeue -u $USER -o "%.10i %.22j %.8T %.10M %.10l %.20R"
