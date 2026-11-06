# FinchExamples

Example problems solved with Finch, a Julia finite element package. The
examples use Finch V2 (`budvinchathura/Finch-V2`, branch `mlir-batched-linalg`)
with its Julia target; the folders in `Reference/` still use the original
[Finch](https://paralab.github.io/Finch/dev/) through `Finch_patched/`.

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
| `Example13_InterfaceNormalFDvsFEM` | The CDI interface normal by finite differences (forward, backward, central) and linear finite elements (element-wise, nodal) on a fixed field |
| `Example17_NavierStokesConservativeDiffusiveInterfaceExplicit2D` | Navier-Stokes + CDI, explicit incremental projection, Finch V2 dev branch: drop advection, static drop, rising bubble (Hysing TC1) |
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
(figure scripts) and `output/` (results; not committed). Finch V2 generates a
Julia project for each run and solves it; that project and its raw results land
in `output/finch/`. In Examples 9 and 10 the scheme is in `scheme.jl`, shared by
`main.jl` and the studies below. Example 9 and References
10, 12 and 14 also have `convergence.sh`, a temporal convergence study
(`convergence.jl`), whose figure is `plot_Figure4_TemporalConvergence.m`. Examples 9
and 10 also have `sharpening_sweep.sh`, a study of how the sharpening term is split
between the unknown and the lagged field (`sharpening_sweep.jl`: five forms, one to
three Picard iterations, six time steps, two meshes; max error at the final time).
Its figures are `plot_Figure5/6` (Example 9) and `plot_Figure6/7` (Example 10), and
`plot_Video2_SharpeningSweep.m` writes one movie per mesh and time step.

Examples 5-8 are explicit (forward Euler). Examples 9-12 are the semi-implicit scheme
(BDF2, two Picard iterations per step), with the sharpening term linearized as in the
derivation, one factor implicit. Finch V2 has forward and backward Euler only, so one
semi-implicit step is a `step_schedule`: phi is solved as a quasi-static equation whose
time derivative is the BDF formula over the previous levels, held as fields, and
`update` entries move the history, the lagged field and the BDF coefficients along.
In References 10 and 11 it is evaluated entirely from
the Picard iterate; References 14 and 15 linearize it by Newton's method about the
Picard iterate, `phi - phi^2 ~ (1 - 2*phis)*phi + phis^2`. References 12 and 13 add the
residual-based VMS stabilization of the derivation to Examples 9 and 10.

Interface normal. Examples 5-12 recover the normal at the nodes by an L2 projection of
the finite element gradient (a continuous field), solved as its own field, one per
component, before each step or Picard iteration. Finch V2 has no lumped mass matrix,
so the projection uses the consistent one; the V1 versions of these examples lumped
it. Compared on 100 x 100 meshes, V1 against V2, the L1 shape error went from 0.0126
to 0.0077 (Example 7), 0.0141 to 0.0076 (Example 8), 0.0130 to 0.0111 (Example 11)
and 0.0142 to 0.0074 (Example 12); Examples 5, 6, 9 and 10 agree with V1 to two to
eight digits, as do the temporal convergence orders of Example 9. References 10-16 still
take it element by element from the gradient of the linear basis, which is constant
per element and jumps between elements; on grid-aligned interfaces that normal grows
a grid-scale zigzag (on the circle of Example 10 it gave a max error of 0.128 with the
semi-implicit scheme, and on Zalesak's disk wiggly slot walls and shed fragments, L1
error 1.2e-2 against 3.4e-3 with the nodal normal). The nodal normal removes it in
every case tested; the Reference examples have not been converted.

Example 7 runs the explicit scheme of Example 6 on Zalesak's slotted disk in
solid-body rotation for one revolution (`elementwise_normal.sh` reruns it with the
element-wise normal, `CDI_NORMAL=elementwise`, results in
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
Example 7 has the same `paper512.sh` for Zalesak's disk: explicit, about 4 hours
(79,616 steps, L1 7.9e-4, phi stays in [0, 1]). These run times and errors are
from the V1 versions; see Run time below for V2.

Example 13 shows why the examples use the nodal normal. On a tanh profile with
eps = h (1D, and a disk in 2D), on the same nodes for every method, it checks
that the linear element-wise gradient is the forward/backward difference of
that element and the nodal gradient is the central difference (exactly in 1D;
in 2D averaged over neighbouring rows with weights 1/6, 2/3, 1/6). Angle error
of the normal on the disk: one-sided and element-wise 20 degrees, central
difference (Mirjalili et al.) 2, nodal 0.07. It takes under a minute.

Example 17 couples the CDI phase field to the incompressible Navier-Stokes
equations (P30 note "NS + CDI, explicit coupling"): forward Euler, an
incremental projection with four linear solves per step (phase field,
density-weighted momentum predictor, variable-coefficient pressure Poisson,
L2 velocity correction) and the consistent momentum flux of Mirjalili and Mani
(2021), U = rho u - (rho1 - rho2) R. It needs Finch V2's development branch
(vector unknowns and fields): run.sh uses ~/Packages/Finch-V2-dev. CDI_CASE
selects the test:
- A, drop advection (rho ratio 10, uniform u = (1, 0), 64^2, 1 min): the
  consistent flux keeps u = (1, 0) to 1e-3 after one period; with
  CDI_THETA23=0 the run blows up at t = 0.94 (run it to CDI_T=0.9 for the plots).
- B, static drop (64^2, 1 min): Laplace jump 3.96 against sigma/R = 4.
- C, rising bubble, Hysing et al. (2009) test case 1 (80 x 160, 17 min): max
  rise velocity 0.2285 (0.2417), min circularity 0.9103 (0.9013), centroid at
  t = 3 1.0618 (1.0813).
  The reference curves and the final shape (TP2D, h = 1/320, files c1g1l7*.txt from
  featflow.de) are in reference_hysing/.
Findings: grad(phi) in the capillary stress must be the L2-projected nodal
gradient (with the element gradient the static drop blows up at t ~ 0.05 for
any dt); the viscous step limit is dx^2/(24 nu) with the consistent mass and
the full stress, tighter than the note's; the singular pressure Poisson of a
closed or periodic box is pinned by a tiny reaction term (CDI_PREG); `dt` in a
quasi-static weak form is not the step size, so it enters as a coefficient.
In case A the velocity error of the consistent run grows from round-off
like e^(11 t), independent of dt: plot_Figure5 shows the pressure carries a
node-to-node (checkerboard) oscillation as large as its smooth part, with the
same pattern as the velocity error and the same growth rate, i.e. a spurious
pressure mode of the equal-order u, p pair that the projection alone does
not suppress. Case C on 120 x 240 (CDI_NEL=120, 1 h, plot_Figure4): 0.2333,
0.9057, 1.0675, all closer to the benchmark.

## Setup, once

1. Install Julia 1.12 (Finch V2 is developed against it), for example with
   juliaup: `juliaup add 1.12`. The run scripts call `julia +1.12`.
2. Clone this repository and Finch V2 (a private repository; ask for access)
   side by side:

       git clone git@github.com:baskargroup/FinchExamples.git
       git clone -b mlir-batched-linalg git@github.com:budvinchathura/Finch-V2.git ~/Packages/Finch-V2

3. Install Finch V2's dependencies (downloads and precompiles, a few minutes):

       cd ~/Packages/Finch-V2
       julia +1.12 --project=. -e 'using Pkg; Pkg.instantiate()'

The run scripts look for Finch V2 in `~/Packages/Finch-V2`; set `FINCH_V2_DIR` to
use another checkout. Only the Julia target is used, so neither IREE nor a C++
toolchain is needed.

## Run

Inside any example folder:

    ./run.sh

Results and a `runlog.txt` with the printed output land in that example's
`output/`. Each run first generates and compiles Finch's code (20-30 s).

Figures need MATLAB: from the example's `matlabscripts/` folder run any
`plot_*.m`; the PDFs are written to `output/` as well.

Run time. Finch V2 re-assembles the nonlinear sharpening term every explicit
step, which makes the explicit 2D examples about five times slower per step than
their V1 versions (about 28 ms per step on 100 x 100): roughly 10 minutes for
Example 6, 45 for Example 7 and 30 for Example 8 at their default meshes, and many
hours for the `paper512.sh` runs. The semi-implicit examples run close to V1's speed
(Example 10 on 100 x 100: 47 s against 36 s). Run one example at a time: two Finch
V2 runs side by side slowed each other down more than tenfold here.

## The Reference examples and `Finch_patched/`

The folders in `Reference/` run on the original Finch, which is not installed as
a Julia package: the repository carries its own copy, upstream
[`paralab/Finch`](https://github.com/paralab/Finch) at commit `defa03c` plus eight
fixes to the FEM boundary-integral code generation. `Finch_patched/PATCHES.md`
lists every change, how it was verified, and what remains unfixed. Finch is MIT
licensed; its license is included. To run a Reference example, install that copy
once (`cd Finch_patched; julia --project=. -e 'using Pkg; Pkg.instantiate()'`,
Julia 1.13) and use its `run.sh`; `FINCH_DIR=/path/to/Finch ./run.sh` runs against
another checkout. Reference 13 is the slow one: about an hour on a laptop (see its
`notes/notes_cost.md`).

## Housekeeping

`./clean.sh` deletes everything git ignores (all results in every `output/`);
`./clean.sh -n` previews what would be removed. Per-folder `.gitignore` files
keep results out of the repository.
