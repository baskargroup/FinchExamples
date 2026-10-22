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
| `Example5_ConservativeDiffusiveInterfaceExplicitTranslating1D` | CDI phase-field model: translating interface |
| `Example6_ConservativeDiffusiveInterfaceExplicitCircle2D` | CDI phase-field model: circle in uniform flow |
| `Example7_ConservativeDiffusiveInterfaceExplicitZalesak2D` | CDI, explicit (as Example 6): Zalesak's slotted disk, one revolution |
| `Example8_ConservativeDiffusiveInterfaceExplicitShearDrop2D` | CDI, explicit (as Example 7): drop in a time-reversed shear (vortex) flow, one period |
| `Example9_ConservativeDiffusiveInterfaceImplicitTranslating1D` | CDI, semi-implicit (BDF2, two Picard iterations): translating interface |
| `Example10_ConservativeDiffusiveInterfaceImplicitCircle2D` | CDI, semi-implicit, as derived: circle in uniform flow |
| `Example11_ConservativeDiffusiveInterfaceImplicitZalesak2D` | CDI, semi-implicit (as Example 10): Zalesak's slotted disk, one revolution |
| `Example12_ConservativeDiffusiveInterfaceImplicitShearDrop2D` | CDI, semi-implicit (as Example 11): drop in a time-reversed shear (vortex) flow, one period |
| `Reference/Example10_ConservativeDiffusiveInterfaceImplicitLagged1D` | as Example 9 with the sharpening term fully lagged |
| `Reference/Example11_ConservativeDiffusiveInterfaceImplicitLagged2D` | as Example 10 with the sharpening term fully lagged |
| `Reference/Example12_ConservativeDiffusiveInterfaceVMS1D` | as Example 9 with VMS stabilization (residual-based, tau from the derivation) |
| `Reference/Example13_ConservativeDiffusiveInterfaceVMS2D` | as Example 10 with VMS stabilization |
| `Reference/Example14_ConservativeDiffusiveInterfaceImplicitNewton1D` | as Example 9 with the sharpening term linearized by Newton's method |
| `Reference/Example15_ConservativeDiffusiveInterfaceImplicitNewton2D` | as Example 10 with the sharpening term linearized by Newton's method |
| `Reference/Example16_ConservativeDiffusiveInterfaceVMSZalesak2D` | as Example 7 with the three VMS stabilization terms (flow and sharpening velocity) |

The folders in `Reference/` (numbered 10-16 there, a numbering of their own that
overlaps the main examples) are earlier variants kept for reference; they are not
being worked on at present. Below, "Reference N" means `Reference/ExampleN_...`.

Each folder contains `main.jl` (the Finch script), `run.sh`, `matlabscripts/`
(figure scripts) and `output/` (results; not committed). Example 9 and References
10, 12 and 14 also have `convergence.sh`, a temporal convergence study
(`convergence.jl`), whose figure is `plot_Figure4_TemporalConvergence.m`. Examples 9
and 10 also have `sharpening_sweep.sh`, a study of how the sharpening term is split
between the unknown and the lagged field (`sharpening_sweep.jl`: five forms, one to
three Picard iterations, six time steps, two meshes; max error at the final time).
Its figures are `plot_Figure5/6` (Example 9) and `plot_Figure6/7` (Example 10), and
`plot_Video2_SharpeningSweep.m` writes one movie per mesh and time step.

Examples 5-8 are explicit (forward Euler). Examples 9-12 are the semi-implicit scheme
(BDF2, two Picard iterations per step), with the sharpening term linearized as in the
derivation, one factor implicit. In References 10 and 11 it is evaluated entirely from
the Picard iterate; References 14 and 15 linearize it by Newton's method about the
Picard iterate, `phi - phi^2 ~ (1 - 2*phis)*phi + phis^2`. References 12 and 13 add the
residual-based VMS stabilization of the derivation to Examples 9 and 10.

Interface normal. Examples 5-12 recover the normal at the nodes by a lumped L2
projection of the finite element gradient (a continuous field). References 10-16 still
take it element by element from the gradient of the linear basis, which is constant
per element and jumps between elements; on grid-aligned interfaces that normal grows
a grid-scale zigzag (on the circle of Example 10 it gave a max error of 0.128 with the
semi-implicit scheme, and on Zalesak's disk wiggly slot walls and shed fragments, L1
error 1.2e-2 against 3.4e-3 with the nodal normal). The nodal normal removes it in
every case tested; the Reference examples have not been converted.

Example 7 runs the explicit scheme of Example 6 on Zalesak's slotted disk in
solid-body rotation for one revolution (about 9 minutes on the 200 x 200 mesh;
`elementwise_normal.sh` reruns it with the element-wise normal, results in
`output/elementwise_normal/`), and
Example 11 runs the semi-implicit scheme on the same problem. Reference 16 adds the
three residual-based stabilization terms of the derivation to Example 7, along the
flow velocity and along the sharpening velocity, all evaluated at the previous level,
with tau computed at the nodes in Julia; the assembly of the stabilization term is
about 50 times more expensive per step than the rest of the scheme, and it still uses
the element-wise normal.

Examples 8 (explicit) and 12 (semi-implicit) solve the drop in shear flow of
Mirjalili et al. (JCP 401, 2020, Section 4.3.1): a drop of radius 0.15 at
(0.5, 0.75) in the vortex of stream function
`(1/pi) sin^2(pi x) sin^2(pi y) cos(pi t/T)`, T = 4. The drop is stretched into a
thin spiral up to T/2 and returns to its initial shape at T, where the L1 (shape)
error is measured. Their Figure 4 shows the most stretched state, t = T/2.
Examples 8 and 12 also have `paper512.sh`, the paper's parameters (eps/dx = 0.7,
gamma/|u|max = 2.5) on a 512 x 512 mesh, results in `output/paper_512/`: explicit
about 4 hours (71,680 steps, L1 5.5e-4), semi-implicit about 2.5 hours (8,192
steps, L1 5.1e-4; phi up to 1.0018 inside the drop while the flow is fast). Set
`CDI_OUTDIR=output/paper_512` before running the MATLAB scripts to plot it.

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
`output/`. The first run of a session takes longer (compilation). Reference 13 is
the slow one: about an hour on a laptop (the stabilized weak form is large; see
its `notes/notes_cost.md`).

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
