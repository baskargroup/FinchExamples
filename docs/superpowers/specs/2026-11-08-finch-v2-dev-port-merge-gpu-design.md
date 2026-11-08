# FinchExamples: port to Finch V2 dev, merge same-equation cases, add a GPU target switch

Date: 2026-11-08. Status: approved in conversation; builder check done, its three
questions settled below as builder calls (marked "Builder call").

## Intent and why

The examples in this repository are the working ground for the CDI (conservative
diffuse interface) and Navier-Stokes + CDI work. Three things are wrong with them
today:

1. Examples 1-16 run on an old Finch V2 branch (`mlir-batched-linalg`, commit
   e5adbc3); Example 17 runs on `dev`. The developer has since folded `dev` into
   `master` and added a GPU path (MLIR/IREE target with `device = :cuda`). Nothing
   in the repository can use it.
2. Several examples solve the same equation with the same scheme and differ only
   by test case (circle, Zalesak, shear drop), each as a separate near-identical
   script. Example 17 shows the better shape: one script, `CDI_CASE` selects the
   test case.
3. The 512^2 runs take hours on a laptop. The GPU path ran a 200^2 Zalesak turn in
   62 s on a Nova A100 (job 17312182, 2026-11-07). Future 3D work needs this.

The outcome wanted: every example on the latest Finch V2, one script per
equation + scheme with a case switch, and the same script able to run on the CPU
(Julia target, the reference) or on a GPU (MLIR target) by a single switch, so
the physics is written once and the two never drift.

## Goals

- All examples run on Finch V2 `origin/dev` at commit 172f4fb (2026-11-06; `master`
  0efb44a contains it). The commit is written down in the README and pinned by
  `FINCH_V2_DIR` defaulting to `~/Packages/Finch-V2`. One checkout with the same
  folder name on both machines: `~/Packages/Finch-V2` (laptop) and
  `/work/mech-ai/suresh/packages/Finch-V2` (Nova), both on branch `dev`. The old
  worktrees (`Finch-V2-dev` locally, `Finch-V2-latest` on Nova) were removed on
  2026-11-08; IREE is installed inside the Nova checkout (`.iree/`).
- Merge only same equation + same scheme + different test case:
  6+7+8 (explicit CDI 2D), 10+11+12 (semi-implicit CDI 2D), 15+16 (CDI + gas 2D).
- Remove Example 13 (interface normal FD vs FEM; it was for a presentation).
- Renumber (table below). Dissolved-gas examples start at 31, Navier-Stokes + CDI
  at 21, so each family has room to grow.
- One file per example with a target switch `FINCH_TARGET` = `julia` (default,
  CPU, the reference) | `mlir_cpu` | `mlir_cuda`. Same physics, same cases, same
  post-processing in every mode.
- Every ported script reproduces its old results before any further change (the
  existing `output/` folders are the regression references).
- Nova job scripts per example (CPU and GPU), modeled on the two that already work
  (`Example17_.../runNova.sh`, `NovaGPU/runNovaGPU.sh`).
- README rewritten for the new numbering, setup (dev branch, IREE on Nova) and the
  target switch.

## Non-goals

- No new physics, no new test cases, no changes to parameters or meshes.
- No GPU version of the Reference/ (Finch V1) examples; they are untouched.
- No change to the MATLAB plotting scripts beyond path/number updates.
- No work inside Finch V2 itself (the Finch folders are read-only for this
  repository); limitations found are reported to the developer, not patched.
- No 3D.

## Anti-goals (what would make this a failure)

- A GPU version whose physics differs from the CPU one, or two copies of a weak
  form that must be edited in step. The switch must only change `target =`,
  the solver settings the MLIR target requires, and the way the run is launched.
- A merged script whose cases silently change results. Each case must reproduce
  the old example's numbers (see Verification).
- Breaking the Zalesak/shear-drop 512^2 results that are in the slide decks.
- Losing the uncommitted Examples 14-16 work (they become 31 and 32).

## Constraints

- Finch V2 MLIR target (as of 172f4fb), confirmed from its code and docs:
  - Solvers: CG or BiCGSTAB only, Jacobi or no preconditioner. Decision: BiCGSTAB
    with Jacobi for every example in MLIR mode. The Julia target keeps whatever
    each example used (LU for 1-4, CG/BiCGSTAB elsewhere); the switch may set the
    solver per target.
  - Fixed dt for the whole run. `dt = T/N` can round to N+1 steps on that target:
    pass `Nsteps` (and `time_interval(T)`), not `dt`.
  - No Julia code during the run: the script emits, then `build.sh` and `run.sh`
    run the compiled program, then Julia post-processes the result files. Both
    targets write `results/<field>_step<k>.bin` + `manifest.json` + `monitors.csv`,
    read by `Finch.read_results_manifest` / `load_snapshot_values`, so the
    post-processing is shared.
  - Schedule is a flat list: a fixed number of Picard iterations = repeated
    entries (as the semi-implicit scripts already do). Quasi-static solves (nodal
    normal projection) are fine. `update` strings may use `ifelse`, `pow`, `old`,
    `dt`, but not `grad`.
  - Coefficients and initial conditions are traced once: arithmetic, trig, tanh,
    sqrt, abs, min, max, ifelse; no Julia `if`, loops or array lookups. (Same rule
    as the Julia target already enforces.)
  - Explicit Euler, implicit Euler, BDF2, theta, RK available; `old(u, k)` and `dt`
    usable in weak forms when the stepper is `STEP_EULER_IMPLICIT` or `STEP_BDF2`
    (explicit/CN steppers rewrite hand-written `old` forms, a known bug: keep
    Example 17's `alphao` via `update`).
  - Vector variables with `g[1]` indexing and per-component Dirichlet; periodic
    meshes; `ds` boundary terms; time-dependent coefficients.
  - `device = :cuda` needs `nvidia-smi` at emit time or `FINCH_CUDA_ARCH` set.
    `setup_iree()` is Linux only: on Nova it is done (Finch-V2-latest/.iree).
    On the Mac only `julia` mode runs; `mlir_*` modes are for Nova (or any Linux
    box with IREE set up).
- Nova: Julia 1.12.7 and the package depot live under
  `/work/mech-ai/suresh/packages` (home quota 10 GB). One Finch clone,
  `/work/mech-ai/suresh/packages/Finch-V2`, branch `dev` at 172f4fb, IREE in its
  `.iree/` (set up by `setup_iree_nova.sh` in that folder, log `setup_iree.log`).
- FinchExamples on Nova is a git clone at
  `/work/mech-ai/suresh/projects/FinchExamples` (done 2026-11-08). Workflow: edit
  and commit on the laptop, push, `git pull` on Nova. Nothing is rsynced any more;
  results (git-ignored `output*/`) stay where they were produced.
- GPU queue: the mech-ai group cap is 15 GPUs; scavenger needs
  `Partition=scavenger QOS=scavenger`; preempted jobs restart from zero (no
  checkpointing on the MLIR target). All our GPU runs are minutes long.
- `scontrol`/`scancel` on the shared cluster are run by the user, not by the
  assistant; `sbatch`, uploads and read-only checks are fine once the user asked
  for the run.
- Examples 14-16 (future 31, 32) and their README text must not be committed
  until the user says so. All other commits are proposed to the user first.
- Finch folders (`~/Packages/Finch*`) are never edited.

## Decided: numbering and contents

| new | folder name | from | cases (`CDI_CASE`) |
|---|---|---|---|
| 1 | Example1_Poisson1D | 1 | one |
| 2 | Example2_Poisson2D | 2 | one |
| 3 | Example3_ConvectionDiffusion1D | 3 | one |
| 4 | Example4_ConvectionDiffusion2D | 4 | one |
| 5 | Example5_ConservativeDiffusiveInterfaceExplicit1D | 5 | translating |
| 6 | Example6_ConservativeDiffusiveInterfaceExplicit2D | 6, 7, 8 | circle, zalesak, sheardrop |
| 7 | Example7_ConservativeDiffusiveInterfaceImplicit1D | 9 | translating; plus convergence.jl and sharpening_sweep.jl drivers |
| 8 | Example8_ConservativeDiffusiveInterfaceImplicit2D | 10, 11, 12 | circle, zalesak, sheardrop; plus sharpening_sweep.jl |
| 21 | Example21_NavierStokesConservativeDiffusiveInterfaceExplicit2D | 17 | A, B, C, D, E (unchanged) |
| 31 | Example31_ConservativeDiffusiveInterfaceScalarTransportExplicit1D | 14 | gascolumn |
| 32 | Example32_ConservativeDiffusiveInterfaceScalarTransportExplicit2D | 15, 16 | bubbles, electrode |

Example 13 is deleted (`git rm`). 1D and 2D stay separate even for the same
equation: the scripts differ in mesh, weak-form syntax and gradient variables,
which is more than a case switch.

## Decided: the target switch

- `FINCH_TARGET` environment variable, read once at the top of `main.jl`:
  `julia` (default) → `JuliaCPUTarget(out_dir = PROJ_DIR, func_name = ...)`;
  `mlir_cpu` → `MLIRTarget(out_dir = PROJ_DIR, func_name = ..., device = :cpu)`;
  `mlir_cuda` → `MLIRTarget(..., device = :cuda)`.
- Everything between `init_finch` and `solve()` is identical in all modes.
- After `solve()`: in `julia` mode `Finch.run_emitted_julia!` as today. In `mlir_*`
  mode the script runs `bash build.sh` and `bash run.sh` in `PROJ_DIR` itself
  (Julia `run(...)`, failing loudly on a non-zero exit), so `./run.sh` is one
  command in every mode. Then the shared post-processing reads the results.
- Solver: the switch sets `solver = :bicgstab, preconditioner = :jacobi` in
  `mlir_*` modes; `julia` mode keeps the example's current choice.
- Time stepping: `time_stepper(...; Nsteps = NSTEPS, ...)` everywhere (both
  targets), never `dt =`.
- `FINCH_CUDA_ARCH` is left to the environment (auto on a GPU node).
- `run.sh` per example stays the entry point (`FINCH_TARGET=mlir_cuda ./run.sh`),
  default `FINCH_V2_DIR=$HOME/Packages/Finch-V2`, uses `julia +1.12`.

## Decided: merged scripts and shared code

- Merged scripts follow Example 17's pattern: `const CASE = get(ENV, "CDI_CASE", <default>)`,
  one `if CASE == ...` block for the case parameters (IC function, velocity
  coefficient, T, BCs, error reference), one shared weak form and schedule.
  Default case: the old lowest-numbered example (circle for 6 and 8, bubbles for 32).
- Output folders: `output/<case>/` for merged scripts (6, 8, 32); single-case
  scripts keep `output/`; Example 21 keeps its `output_case<X>` folders untouched
  (the slide-deck figure scripts read them). `CDI_OUTDIR` still overrides.
- `common.jl` at the repository root, included by the CDI scripts: snapshot
  reading into time series CSVs, lumped-weight L1 error, mass change from
  monitors.csv, the nodal-normal expression string, and the target switch
  helper itself (`finch_target(name, proj_dir)` returning the target and
  solver keywords, plus `finch_run!(ctx, vars...)` that dispatches on the mode).
  Examples 1-4 may include it for the switch only.
- Dev-branch simplifications applied during the port (results must not change
  beyond round-off): `old(phi, k)` + `dt` replace the `a0, a1, a2, phio, phioo`
  fields of the semi-implicit scripts (`STEP_BDF2` declared for the loop;
  BDF1 first step as before); one vector variable `g` with `g[1]`, `g[2]`
  replaces `gx, gy` + `ex`/`ey` in 2D; `dt` replaces Example 17's `dtc`/`idt`.
  Picard stays as repeated schedule entries.
  Builder call (from the builder check): `old(phi)` inside phi's own equation is
  the current Picard iterate, not phi^n, so the previous level `phio` stays an
  `update`'d field; only `phioo` (-> `old(phi, 2)`) and the BDF coefficient fields
  `a0, a1, a2` (-> `ifelse(t < 1.5*dt, ...)` or a t-dependent coefficient) are
  replaced. The `scheme.jl` `iters` loop and the sweep/convergence drivers keep
  their structure.
  Builder call: Example 21 keeps its `dtc`/`idt` coefficients. A bare `dt` would
  turn its quasi-static solves into time steps under `STEP_EULER_EXPLICIT`.
  Builder call: the case block of a merged script may switch term strings
  (the circle's non-conservative advection vs the conservative split form),
  boundary calls, schedule entries (the electrode's `psi` solve) and monitors,
  exactly as Example 17's cases switch flags and BCs. No case is rewritten to a
  common form; every case reproduces the old example's numbers.
- `paper512.sh`, `elementwise_normal.sh`, `convergence.sh`, `sharpening_sweep.sh`
  are kept and point at the merged scripts with the right `CDI_CASE`.
- MATLAB scripts move with their example; only paths and numbers change.

## Decided: Nova

- Per example: `runNova.sh` (CPU, Julia target, pattern of the Example 17 one) and
  `runNovaGPU.sh` (`FINCH_TARGET=mlir_cuda`, A100 header, cmake/ninja/cuda
  modules, pattern of `NovaGPU/runNovaGPU.sh`). Both bash, both read `CDI_CASE`
  from the environment, both write `run_<jobid>.txt` beside the script.
- The repository is a git clone on Nova at
  `/work/mech-ai/suresh/projects/FinchExamples` (done); `git pull` there, commits
  and pushes only from the laptop. Job scripts are committed so they arrive by pull.
- `NovaGPU/` (the throwaway Finch-examples GPU test) is removed once Example 6's
  GPU run exists.

## Decided: order of work and verification

1. Port each existing script to dev as is (no merge yet), run it on the laptop
   in `julia` mode, compare the printed summary numbers with the old
   `output/runlog.txt` (max error, L1, mass change, min/max phi, run-specific
   checks such as Example 17's impact time). Tolerance: 1e-8 relative where the
   arithmetic is unchanged; solver-tolerance level (1e-6) where LU became a
   Krylov solve. Record both numbers in the README of the example.
   Long runs (the 300^2 circle, the 512^2 cases) are compared at reduced mesh
   on the laptop and at full size on Nova.
2. Merge 6+7+8, 10+11+12, 15+16; re-run every case; same comparison.
3. Delete Example 13, renumber folders, update README, paper512.sh and friends,
   MATLAB paths. Old numbers in slide decks and in the assistant's memory notes
   are not touched; a mapping table goes in the README.
4. GPU: Example 6 zalesak first (Nova, `mlir_cuda`), compared with the CPU
   result at the same mesh, then 512^2 against `paper_512`. Then 1-5, 7, 8, 21,
   31, 32, each compared with its CPU run. `mlir_cpu` is run once per example
   on Nova as the cheaper cross-check.
5. Commits proposed to the user at the end of each step; 31/32 excluded until
   released.

## Left to the builder

- Exact layout of `common.jl` and the helper names, as long as the scripts read
  like Example 17.
- How the 1D scripts express the gradient variable (scalar `gx` stays fine).
- Whether Example 21's `dump_snapshots.jl` survives or folds into `common.jl`.
- Slurm resource numbers per example (time, memory), following the measured runs.
- Any further duplicate code worth sharing beyond the listed blocks.
