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
| `Example7_ConservativeDiffusiveInterfaceImplicit1D` | CDI, semi-implicit (BDF2, two Picard iterations): translating interface |
| `Example8_ConservativeDiffusiveInterfaceImplicit2D` | CDI, semi-implicit, as derived: circle in uniform flow |
| `Example9_ConservativeDiffusiveInterfaceImplicitLagged1D` | as Example 7 with the sharpening term fully lagged |
| `Example10_ConservativeDiffusiveInterfaceImplicitLagged2D` | as Example 8 with the sharpening term fully lagged |
| `Example11_ConservativeDiffusiveInterfaceVMS1D` | as Example 7 with VMS stabilization (residual-based, tau from the derivation) |

Each folder contains `main.jl` (the Finch script), `run.sh`, `matlabscripts/`
(figure scripts) and `output/` (results; not committed). Examples 7, 9 and 11 also have
`convergence.sh`, a temporal convergence study (`convergence.jl`), whose figure is
`plot_Figure4_TemporalConvergence.m`.

Examples 7-10 are the semi-implicit scheme (BDF2, two Picard iterations per step).
In 7 and 8 the sharpening term is linearized as in the derivation, with one factor
implicit; in 9 and 10 it is evaluated entirely from the Picard iterate. The two
agree in 1D, but on the 300 x 300 mesh the derived form grows a grid-scale
wrinkle along the leading edge (`plot_Figure5_FrontZoom.m` in Examples 8 and 10
shows it), which the fully lagged form does not have.

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
