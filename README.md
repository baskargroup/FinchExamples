# FinchExamples

Example problems solved with [Finch](https://paralab.github.io/Finch/dev/), a
Julia finite element package.

## Examples

| Folder | Problem |
|---|---|
| `Example1_Poisson1D` | Poisson equation, 1D |
| `Example2_Poisson2D` | Poisson equation, 2D |
| `Example3_ConvectionDiffusion1D` | Convection-diffusion, 1D |
| `Example4_ConvectionDiffusion2D` | Convection-diffusion, 2D |
| `Example5_ConservativeDiffusiveInterface1D` | CDI phase-field model: translating interface |
| `Example6_ConservativeDiffusiveInterface2D` | CDI phase-field model: circle in uniform flow |

Each folder contains `main.jl` (the Finch script), `run.sh`, `matlabscripts/`
(figure scripts) and `output/` (results; not committed).

## Setup, once

1. Install Julia. These examples were run with Julia 1.13.
2. Clone the repository:

       git clone git@github.com:baskargroup/FinchExamples.git

3. Install Finch's dependencies into the copy of Finch that ships with the
   repository (downloads and precompiles, a few minutes):

       cd FinchExamples/Finch_patched
       julia --project=. -e 'using Pkg; Pkg.instantiate()'

   A message that the optional PETSc package is not available is
   informational, not an error; Finch's built-in solvers do not need it.

## Run

Inside any example folder:

    ./run.sh

Results and a `runlog.txt` with the printed output land in that example's
`output/`. The first run of a session takes longer (compilation).

Figures need MATLAB: from the example's `matlabscripts/` folder run any
`plot_*.m`; the PDFs are written to `output/` as well.

## About `Finch_patched/`

Finch is not installed as a Julia package. The repository carries its own copy,
upstream [`paralab/Finch`](https://github.com/paralab/Finch) at commit
`defa03c` plus eight fixes to the FEM boundary-integral code generation that
Example 6 needs for its outflow boundary condition. `Finch_patched/PATCHES.md`
lists every change, how it was verified, and what remains unfixed. Finch is MIT
licensed; its license is included.

To run against a different Finch checkout instead:

    FINCH_DIR=/path/to/Finch ./run.sh

## Housekeeping

`./clean.sh` deletes everything git ignores (all results in every `output/`);
`./clean.sh -n` previews what would be removed. Per-folder `.gitignore` files
keep results out of the repository.
