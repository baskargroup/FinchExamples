# FinchExamples

Example problems solved with Finch, a Julia finite element package. The
examples use Finch V2 (`budvinchathura/Finch-V2`, branch `dev`, commit
`172f4fb`, 2026-11-06). Every example is one script that runs on the CPU with
Finch's Julia target (the reference) or, with `FINCH_TARGET=mlir_cuda`, on an
NVIDIA GPU through Finch's MLIR/IREE target (see Run). The folders in
`Reference/` still use the original [Finch](https://paralab.github.io/Finch/dev/)
through `Finch_patched/`.

## Examples

| Folder | Problem | Cases (`CDI_CASE`) |
|---|---|---|
| `Example1_Poisson1D` | Poisson equation, 1D | |
| `Example2_Poisson2D` | Laplace equation, 2D, quadratic elements | |
| `Example3_ConvectionDiffusion1D` | Convection-diffusion, 1D | |
| `Example4_ConvectionDiffusion2D` | Convection-diffusion, 2D, quadratic elements | |
| `Example5_ConservativeDiffusiveInterfaceExplicit1D` | CDI phase-field model, explicit (forward Euler): translating interface | |
| `Example6_ConservativeDiffusiveInterfaceExplicit2D` | CDI, explicit, 2D | `circle` in uniform flow, `zalesak` (slotted disk, one revolution), `sheardrop` (drop in a time-reversed vortex, one period) |
| `Example7_ConservativeDiffusiveInterfaceImplicit1D` | CDI, semi-implicit (BDF2, two Picard iterations): translating interface; temporal convergence and sharpening-term studies | |
| `Example8_ConservativeDiffusiveInterfaceImplicit2D` | CDI, semi-implicit, 2D; sharpening-term study | `circle`, `zalesak`, `sheardrop` |
| `Example21_NavierStokesConservativeDiffusiveInterfaceExplicit2D` | Navier-Stokes + CDI, explicit incremental projection | `A` drop advection, `B` static drop, `C` rising bubble (Hysing TC1), `D` dam break (Martin and Moyce), `E` falling drop |
| `Example31_ConservativeDiffusiveInterfaceScalarTransportExplicit1D` | CDI with dissolved-gas transport, explicit: a gas column grows in a supersaturated liquid (Qiu and Mani) | |
| `Example32_ConservativeDiffusiveInterfaceScalarTransportExplicit2D` | CDI with dissolved-gas transport, explicit, 2D | `bubbles` (seeded bubbles in a supersaturated liquid, u = 0), `electrode` (bubbles on a gas-producing electrode, expansion flow) |
| `Reference/Example10_ConservativeDiffusiveInterfaceImplicitLagged1D` | as Example 7 with the sharpening term fully lagged | |
| `Reference/Example11_ConservativeDiffusiveInterfaceImplicitLagged2D` | as Example 8 (circle) with the sharpening term fully lagged | |
| `Reference/Example12_ConservativeDiffusiveInterfaceVMS1D` | as Example 7 with VMS stabilization (residual-based, tau from the derivation) | |
| `Reference/Example13_ConservativeDiffusiveInterfaceVMS2D` | as Example 8 (circle) with VMS stabilization | |
| `Reference/Example14_ConservativeDiffusiveInterfaceImplicitNewton1D` | as Example 7 with the sharpening term linearized by Newton's method | |
| `Reference/Example15_ConservativeDiffusiveInterfaceImplicitNewton2D` | as Example 8 (circle) with the sharpening term linearized by Newton's method | |
| `Reference/Example16_ConservativeDiffusiveInterfaceVMSZalesak2D` | as Example 6 (zalesak) with the three VMS stabilization terms (flow and sharpening velocity) | |

The numbers leave room in each family: 1-8 are the scalar and CDI transport
examples, 21- Navier-Stokes + CDI, 31- CDI with dissolved gas. Examples that
solve the same equation with the same scheme on different test problems are
one script with a case switch, as Example 21 always was. Before 2026-11-08 the
cases were separate examples:

| old | now |
|---|---|
| 5 | 5 |
| 6, 7, 8 | 6 (`circle`, `zalesak`, `sheardrop`) |
| 9 | 7 |
| 10, 11, 12 | 8 (`circle`, `zalesak`, `sheardrop`) |
| 13 (interface normal FD vs FEM) | removed; it was made for a presentation and its figures live in the slide decks |
| 14 | 31 |
| 15, 16 | 32 (`bubbles`, `electrode`) |
| 17 | 21 |

The folders in `Reference/` (numbered 10-16 there, a numbering of their own that
overlaps the main examples) are earlier variants kept for reference; they are not
being worked on at present. Below, "Reference N" means `Reference/ExampleN_...`.

Each folder contains `main.jl` (the Finch script), `run.sh`, `matlabscripts/`
(figure scripts) and `output/` (results; not committed; one subfolder per case
in the case-switched examples, `output/<case>/`, and `output_case<X>/` in
Example 21). Finch V2 generates a project for each run and solves it; that
project and its raw results land in the `finch/` subfolder of the output
folder. `common.jl` at the top level holds what every script shares: the
`FINCH_TARGET` switch, the launch of the generated project, reading the
snapshots back, the lumped L1 error, the mass monitor and the interface-normal
expression. In Examples 7 and 8 the scheme is in `scheme.jl`, shared by
`main.jl` and the studies below. Example 7 and References 10, 12 and 14 also
have `convergence.sh`, a temporal convergence study (`convergence.jl`), whose
figure is `plot_Figure4_TemporalConvergence.m`. Examples 7 and 8 also have
`sharpening_sweep.sh`, a study of how the sharpening term is split between the
unknown and the lagged field (`sharpening_sweep.jl`: five forms, one to three
Picard iterations, six time steps, two meshes; max error at the final time, on
the circle in Example 8). Its figures are `plot_Figure5/6` (Example 7) and
`plot_Figure6/7` (Example 8), and `plot_Video2_SharpeningSweep.m` writes one
movie per mesh and time step.

Examples 5 and 6 are explicit (forward Euler). Examples 7 and 8 are the
semi-implicit scheme (BDF2, two Picard iterations per step), with the
sharpening term linearized as in the derivation, one factor implicit. One
semi-implicit step is a `step_schedule`: phi is solved as a quasi-static
equation whose time derivative is the BDF formula over the previous levels
(the previous level a field moved by an `update`, the one before it Finch's
own history `old(phi, 2)`, the BDF coefficients coefficients of t that switch
from BDF1 to BDF2 after the first step), and `update` entries move the lagged
field along. In References 10 and 11 the sharpening term is evaluated entirely
from the Picard iterate; References 14 and 15 linearize it by Newton's method
about the Picard iterate, `phi - phi^2 ~ (1 - 2*phis)*phi + phis^2`. References
12 and 13 add the residual-based VMS stabilization of the derivation to
Examples 7 and 8.

Interface normal. Examples 5-8 recover the normal at the nodes by an L2
projection of the finite element gradient (a continuous field), solved as its
own vector field before each step or Picard iteration. Finch V2 has no lumped
mass matrix, so the projection uses the consistent one; the V1 versions of
these examples lumped it. Compared on 100 x 100 meshes, V1 against V2, the L1
shape error went from 0.0126 to 0.0077 (Example 6 zalesak), 0.0141 to 0.0076
(Example 6 sheardrop), 0.0130 to 0.0111 (Example 8 zalesak) and 0.0142 to
0.0074 (Example 8 sheardrop); the 1D cases and the circles agree with V1 to two
to eight digits, as do the temporal convergence orders of Example 7. References
10-16 still take it element by element from the gradient of the linear basis,
which is constant per element and jumps between elements; on grid-aligned
interfaces that normal grows a grid-scale zigzag (on the circle of Example 8 it
gave a max error of 0.128 with the semi-implicit scheme, and on Zalesak's disk
wiggly slot walls and shed fragments, L1 error 1.2e-2 against 3.4e-3 with the
nodal normal). The nodal normal removes it in every case tested
(`elementwise_normal.sh` in Example 6 reruns Zalesak's disk with the
element-wise normal, `CDI_NORMAL=elementwise`); the Reference examples have not
been converted. Reference 16 adds the three residual-based stabilization terms
of the derivation to the explicit Zalesak case, along the flow velocity and
along the sharpening velocity, all evaluated at the previous level, with tau
computed at the nodes in Julia; the assembly of the stabilization term is about
50 times more expensive per step than the rest of the scheme, and it still uses
the element-wise normal.

The `sheardrop` cases of Examples 6 (explicit) and 8 (semi-implicit) solve the
drop in shear flow of Mirjalili et al. (JCP 401, 2020, Section 4.3.1): a drop
of radius 0.15 at (0.5, 0.75) in the vortex of stream function
`(1/pi) sin^2(pi x) sin^2(pi y) cos(pi t/T)`, T = 4. The drop is stretched into
a thin spiral up to T/2 and returns to its initial shape at T, where the L1
(shape) error is measured. Their Figure 4 shows the most stretched state,
t = T/2. `paper512_sheardrop.sh` (Example 6) and `paper512.sh` (Example 8) run
the paper's parameters (eps/dx = 0.7, gamma/|u|max = 2.5) on a 512 x 512 mesh,
results in `output/paper_512/`: explicit about 4 hours (71,680 steps, L1
5.5e-4), semi-implicit about 2.5 hours (8,192 steps, L1 5.1e-4; phi up to
1.0018 inside the drop while the flow is fast). Set
`CDI_OUTDIR=output/paper_512` before running the MATLAB scripts to plot it.
`paper512_zalesak.sh` (Example 6) does the same for Zalesak's disk: explicit,
about 4 hours (79,616 steps, L1 7.9e-4, phi stays in [0, 1]).

Example 31 adds a second unknown, the dissolved gas concentration c, to the
explicit phase field of Example 5. Gas crossing the interface at the rate
J = alpha/tau max(0, c - csat(1 - alpha)) expands the column. One step is a
`step_schedule`: update J; solve psi_xx = theta15 J (psi = 0 at both ends) for
the velocity u = psi_x; project grad(alpha) to the nodes for the normal; step
c; step alpha. A column of width 1/4 in liquid at 3 csat grows to
[0.294, 0.706] by t = 1 (200 elements, about 25 s); the gas volume balance
closes to 4e-7. Two properties of the model show up: trace gas in the liquid
grows like exp(2 theta15 t/tau) because J is proportional to alpha, so small gas
pockets (alpha up to ~0.05) form ahead of the interfaces, and with eps = dx the
growth rate depends on the mesh (the transfer acts only inside the interface, so
it scales with eps; 400 elements give [0.331, 0.669]). The pockets are removed
by a nucleation threshold, J = 0 where alpha <= 0.01 (CDI_AMIN, default 0.01;
0 is the note's model): alpha ahead of the column stays below 1e-9 and the
column edges move by 6e-4, to [0.2946, 0.7054].

Example 32 `bubbles` is the 2D case of the same note: five bubbles (radii 0.08,
0.06, 0.04 isolated, and a pair of 0.06 with a 0.06 gap) in [0, 1]^2, alpha = 0
and c = 3 csat on the boundary, u = 0, tau = 0.5, D = 0.01. The normal is the
nodal one, and J has the nucleation threshold of Example 31. On 100 x 100 to
t = 1.5 (110 s) the volume balance closes to 5e-6 and alpha stays in [0, 1];
the radii roughly double, and the pair merges at t ~ 1 and then joins bubble
B. With u = 0 no liquid is displaced, so the dissolved gas of the liquid the
bubbles take over is pushed outward and concentrated: c rises above csuper (to
6.3) around the bubbles instead of being depleted, and the growth speeds up.
t = 1.5 is about as far as this case goes: past it gas is created inside
bubbles already at alpha = 1, alpha exceeds 1 and the explicit step blows up
near t = 2.2. The default 200 x 200 run takes about 40-50 minutes.

Example 32 `electrode` puts nucleation sites (half discs) on a gas-producing
electrode: the bottom wall injects dissolved gas at a constant rate, the
liquid starts saturated, the sides are periodic and the top is saturated
liquid. The velocity is the expansion flow u = grad(psi), lap(psi) = theta15 J,
so the growing bubbles push the liquid out through the top. On 100 x 100
(4.5 min) the gas volume grows from 0.022 to 0.149 by t = 2, balanced to 2e-8.
Without buoyancy or surface tension a wall bubble grows sideways as fast as
upward, so the bubbles merge into a gas film along the wall by t ~ 1.3. The
default 200 x 200 run takes about 35 minutes.

Example 21 couples the CDI phase field to the incompressible Navier-Stokes
equations (P30 note "NS + CDI, explicit coupling"): forward Euler, an
incremental projection with four linear solves per step (phase field,
density-weighted momentum predictor, variable-coefficient pressure Poisson,
L2 velocity correction) and the consistent momentum flux of Mirjalili and Mani
(2021), U = rho u - (rho1 - rho2) R. CDI_CASE selects the test:
- A, drop advection (rho ratio 10, uniform u = (1, 0), 64^2, 1 min): the
  consistent flux keeps u = (1, 0) to 1e-3 after one period; with
  CDI_THETA23=0 the run blows up at t = 0.94 (run it to CDI_T=0.9 for the plots).
- B, static drop (64^2, 1 min): Laplace jump 3.96 against sigma/R = 4.
- C, rising bubble, Hysing et al. (2009) test case 1 (80 x 160, 17 min): max
  rise velocity 0.2285 (0.2417), min circularity 0.9103 (0.9013), centroid at
  t = 3 1.0618 (1.0813).
  The reference curves and final shapes (TP2D levels 4-7, h = 1/40 to 1/320,
  files c1g1l<n>*.txt from featflow.de) are in reference_hysing/; plot_Figure6
  compares at equal h.
- D, dam break, Martin and Moyce (1952) column a x 2a in a 4a x 3a box, free
  slip walls, non-dimensional with the box width and g = 1, water/air density
  and surface tension (rho 1/0.001, sigma 2.2e-5) but the water viscosity ten
  times the physical one (mu 1e-3/1.8e-6, Re ~ 1500), 128 x 96, 25 min: the
  front reaches the far wall at T = t sqrt(2g/a) = 3.02 (experiment ~3) and
  runs 5-10% ahead of the measured points before that. front.csv holds the
  front position and column height; the experimental front (digitised as in
  the Lethe dam-break example) is in reference_martin_moyce/. With the
  physical viscosity (CDI_MU1=1e-4) a mode grows on the air side of the column
  face and the run blows up at t ~ 0.2 for any dt and for density ratios 10 to
  1000 (the same equal-order pressure mode as in case A, amplified by dt/rho in
  the air); dump_snapshots.jl writes the snapshots of such a run for plotting.
- E, falling drop: a water drop (R = 0.1) released at rest in air at (0.5, 1.5)
  above a pool of depth 0.5 in [0,1]x[0,2], free slip walls, the fluids of
  case D but water viscosity 5e-3, Bond number 2, 96 x 192, 37 min: drop.csv
  holds the centroid and mean vertical velocity of the water above the pool;
  they follow free fall (y_0 - g t^2/2, -g t) to contact at t = 1.33, then a
  qualitative splash (entrained bubbles, crater, rebound) to t = 2.5. With the
  dam-break viscosity 1e-3 (CDI_MU1=1e-3) the fall is identical but the impact
  grows the air-side mode and the run blows up at t ~ 1.46.
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
   side by side, Finch V2 on its `dev` branch:

       git clone git@github.com:baskargroup/FinchExamples.git
       git clone -b dev git@github.com:budvinchathura/Finch-V2.git ~/Packages/Finch-V2

   The examples were last verified at Finch V2 commit `172f4fb` (2026-11-06);
   `git -C ~/Packages/Finch-V2 checkout 172f4fb` pins it.
3. Install Finch V2's dependencies (downloads and precompiles, a few minutes):

       cd ~/Packages/Finch-V2
       julia +1.12 --project=. -e 'using Pkg; Pkg.instantiate()'

The run scripts look for Finch V2 in `~/Packages/Finch-V2`; set `FINCH_V2_DIR` to
use another checkout. The Julia target needs nothing else. The MLIR target
(`FINCH_TARGET=mlir_cpu` or `mlir_cuda`) needs IREE, installed once on a Linux
machine with internet, cmake, ninja and a C++ compiler (and the CUDA toolkit
for the GPU) by `julia --project=. -e 'using Finch; setup_iree(device = :cuda)'`
in the Finch checkout; Finch's README has the details. On Nova that is done
(`NOVA_SSH_SETUP_MAC.md` describes the Nova set-up).

## Run

Inside any example folder:

    ./run.sh                              CPU, Julia target: the reference
    CDI_CASE=zalesak ./run.sh             a case of a case-switched example
    FINCH_TARGET=mlir_cuda ./run.sh       the same run on the GPU (Linux + IREE)
    FINCH_TARGET=mlir_cpu ./run.sh        the MLIR target on the CPU

Results and a `runlog.txt` with the printed output land in that example's
`output/` (`output/<case>/`, or `output_case<X>/` in Example 21). The physics
(mesh, variables, weak forms, schedule) is the same in every mode; the switch
changes only the code target, the solver (the MLIR target has CG and BiCGSTAB
with Jacobi, so it uses BiCGSTAB everywhere) and how the generated project is
launched. In `mlir_*` modes the script emits the project, runs its `build.sh`
(iree-compile and cmake) and `run.sh`, and post-processes the results in
Julia afterwards. Each run first generates (and in Julia mode compiles) Finch's
code, 20-30 s.

Figures need MATLAB: from the example's `matlabscripts/` folder run any
`plot_*.m`; the PDFs are written to the output folder as well.

Run time (laptop, Julia target). Finch V2 re-assembles the nonlinear sharpening
term every explicit step: Example 6 `circle` about 9 minutes, `zalesak` about
45 and `sheardrop` about 30 at their default meshes, and hours for the
`paper512*.sh` runs. The semi-implicit cases run in minutes. Run one example at
a time: two Finch V2 runs side by side slow each other down more than tenfold.

## Reference numbers

`regression/` holds the summary lines each example printed when it was last
verified, and `tools/compare_runlog.py REF NEW` checks a new `runlog.txt`
against them (relative tolerance 1e-8 on the scale of each line, absolute
1e-12). After any change to a script or to Finch:

    ./run.sh && python3 ../tools/compare_runlog.py ../regression/<reference>.txt output/runlog.txt

The references were taken on the Julia target; the MLIR target solves with
BiCGSTAB instead of LU or CG, so its numbers agree to the solver tolerance,
not to 1e-8. `bash tools/test_compare_runlog.sh` and
`julia +1.12 --project=$HOME/Packages/Finch-V2 tools/test_common.jl` test the
tools themselves.

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
