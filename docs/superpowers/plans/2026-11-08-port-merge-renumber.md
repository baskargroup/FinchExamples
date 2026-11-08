# Port to Finch V2 dev, merge same-equation cases, renumber: Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Every example runs on Finch V2 dev 172f4fb through one `FINCH_TARGET` switch, the three same-equation groups are merged into case-switched scripts, folders are renumbered, and every ported case reproduces its old numbers.

**Architecture:** One `common.jl` at the repository root holds the target switch and the post-processing blocks that every CDI script copies today; each example's `main.jl` keeps Example 17's shape (env-selected case, case block, one weak form + schedule, shared post-processing). Regression is a small Python comparer run on `output/runlog.txt` against committed reference numbers taken from the old runs. GPU verification and Nova job scripts are a second plan, after this one (the switch's `mlir_*` modes are implemented here but only exercised in `julia` mode).

**Tech Stack:** Julia 1.12 (`julia +1.12`), Finch V2 dev at `~/Packages/Finch-V2` (branch `dev`, 172f4fb), zsh run scripts, Python 3 for the comparer, MATLAB plot scripts (paths only).

**Spec:** `docs/superpowers/specs/2026-11-08-finch-v2-dev-port-merge-gpu-design.md`

## Global Constraints

- Finch checkout: `FINCH_V2_DIR` default `$HOME/Packages/Finch-V2`, commit 172f4fb; the Finch folders are never edited.
- Time stepping declared with `Nsteps = NSTEPS` (and `time_interval(T)`), never `dt =`.
- `FINCH_TARGET` values: `julia` (default), `mlir_cpu`, `mlir_cuda`. `julia` mode keeps each example's current solver (LU for Examples 1-4); `mlir_*` modes set `solver = :bicgstab, preconditioner = :jacobi`.
- Case variable `CDI_CASE`; output of merged scripts in `output/<case>/`; `CDI_OUTDIR` overrides; Example 21 keeps `output_case<X>`.
- Reproduction tolerance: relative 1e-8 on every printed summary number (same arithmetic, same solver); a failure is reported to the user, never patched into the reference.
- Run one Finch job at a time on the laptop (README: side-by-side runs slow each other tenfold).
- Examples 14, 15, 16 (future 31, 32) and README text about them are not committed until the user says so. Other commits are made as the plan says; `git push` only in the last task.
- No new physics, parameters or meshes.

## Review Focus

1. `FINCH_TARGET=foo` must stop with a message listing the three values before any Finch call (test in Task 2).
2. `CDI_CASE=foo` on a merged script must stop with a message listing the cases (tests in Tasks 6, 8, 10).
3. Running two cases of a merged script must not overwrite each other: `output/circle/` and `output/zalesak/` both present afterwards (test in Task 6).
4. The printed `dt` after the `Nsteps` change must equal the old printed `dt` to 1e-12, otherwise every later number drifts (checked by the comparer in Tasks 4-11, `dt` is one of the compared numbers).
5. `mlir_*` on a machine without IREE must fail at `build.sh` with its own message, not with a Julia stack trace from the post-processing reading missing results (test in Task 2, mocked: `FINCH_TARGET=mlir_cpu` on Example 1 must exit non-zero and print `[build.sh]`).

---

### Task 1: Regression comparer and reference numbers

**Files:**
- Create: `tools/compare_runlog.py`
- Create: `regression/<old folder name>.txt` for Examples 1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 15, 17 (case A) — the summary lines copied from the existing `output/runlog.txt` (Example 2 prints no numbers: skip; Example 16 has no log: Task 9 makes it)
- Create: `tools/test_compare_runlog.sh`

**Interfaces:**
- Produces: `python3 tools/compare_runlog.py REF NEW [--rtol 1e-8]` exit 0 when every number in REF is found in NEW within tolerance, exit 1 with a table of mismatches otherwise. A "number" is every float in a line of the form `name = v`, `name: a -> b`, `name: a -> b, relative change c`, `dt = v, steps = n ...`, `bubble X: ... r1 -> r2`, `column: [a, b] -> [c, d]`; lines are matched by their text before the first `=` or `:`; lines of REF absent from NEW are mismatches; lines starting with `run time` are ignored.

- [ ] **Step 1: Write the test** `tools/test_compare_runlog.sh`: creates two temp files with the Example 7 summary lines (from this plan's Task 6 table), runs the comparer on identical files (expect exit 0), on a copy with `L1 error = 0.003385368754585961` changed to `0.0033853687` (expect exit 1 and the word `L1` in the output), on a copy with `run time = 1 s` added (expect exit 0), and on a copy missing the `max phi` line (expect exit 1).
- [ ] **Step 2: Run it**: `bash tools/test_compare_runlog.sh` — expected: fails with "No such file" for the comparer.
- [ ] **Step 3: Implement `tools/compare_runlog.py`** (stdlib only: `re`, `sys`, `argparse`). Parse both files into `{key: [floats]}` where key = the line text before the first `=` or `:` with trailing spaces stripped; compare lists pairwise with `abs(a-b) <= rtol*max(|a|,|b|)` and a 1e-300 absolute floor for zeros; print `OK (N numbers)` or one line per mismatch `key  index  ref  new  rel`.
- [ ] **Step 4: Run the test**: `bash tools/test_compare_runlog.sh` — expected: `ALL PASS`.
- [ ] **Step 5: Extract the references**: for each example folder listed above, `grep -v -E '^\[|^ *@|^\s*$|Precompil|✓|Finch\]|^ [├└─ ]|^─|Finch has completed' Example*/output/runlog.txt` keeping only the summary lines (the block from the `dt = ...` or first `max error` line to `N snapshots written`), into `regression/<folder>.txt`. For Example 17 use `output_caseA/runlog.txt` into `regression/Example17_..._caseA.txt`.
- [ ] **Step 6: Check each reference against its own runlog**: `for f in regression/*.txt; do python3 tools/compare_runlog.py $f <matching runlog>; done` — expected: every line `OK`.
- [ ] **Step 7: Commit**: `git add tools regression && git commit -m "Regression comparer and reference numbers from the old runs"` (do not add Examples 14-16 references yet: keep `regression/Example14*.txt`, `Example15*.txt` untracked until released).

### Task 2: `common.jl` with the target switch, run dispatch and shared post-processing

**Files:**
- Create: `common.jl`
- Create: `tools/test_common.jl`
- Modify: `Example1_Poisson1D/main.jl`, `Example1_Poisson1D/run.sh` (first user of the switch)

**Interfaces:**
- Produces, in `common.jl` (included with `include(joinpath(@__DIR__, "..", "common.jl"))`):
  - `const FINCH_TARGET = get(ENV, "FINCH_TARGET", "julia")`; any other value than `julia|mlir_cpu|mlir_cuda` → `error("FINCH_TARGET must be julia, mlir_cpu or mlir_cuda, got ...")` at include time.
  - `finch_target(func_name::String, proj_dir::String) -> Finch.AbstractTarget` (`JuliaCPUTarget(out_dir = proj_dir, func_name = func_name)` or `MLIRTarget(out_dir = proj_dir, func_name = func_name, device = :cpu | :cuda)`).
  - `finch_solver_kwargs(julia_solver::Symbol; tol = 1e-9, max_iter = 500) -> NamedTuple` → `(solver = julia_solver, solver_tol = tol, solver_max_iter = max_iter)` in julia mode, `(solver = :bicgstab, preconditioner = :jacobi, solver_tol = tol, solver_max_iter = max_iter)` otherwise; scripts splat it into `init_finch`.
  - `finch_run!(ctx, vars...)`: julia mode → `Finch.run_emitted_julia!(ctx.target, vars...)`; mlir modes → `run(Cmd(\`bash build.sh\`; dir = proj_dir))` then `run(Cmd(\`bash run.sh\`; dir = proj_dir))` where `proj_dir = ctx.target.out_dir`; a non-zero exit raises (Julia `run` does).
  - `read_snapshots(proj_dir::String) -> (times::Vector{Float64}, fields::Dict{String, Matrix{Float64}})`, matrix = nodes × snapshots, from `Finch.read_results_manifest` / `load_snapshot_values` on `joinpath(proj_dir, "results")`.
  - `final_values(fields, name, N) -> Vector` = last column; vector fields are component-blocked (`v[1:N]`, `v[N+1:2N]`), helper `final_vector(fields, name, N) -> (c1, c2)`.
  - `write_timeseries(out_dir, name, coords::Matrix, A::Matrix)` → `joinpath(out_dir, name * "_timeseries.csv")` with `[coords A]`, and `write_snapshot_times(out_dir, times)`.
  - `lumped_weights(x, y, dx) -> Vector` (dx² inside, halved on edges, quartered at corners: exactly the loop of Example 7 lines 190-194).
  - `mass_change(results_dir, col = 3) -> (m0, m1, rel)` from `monitors.csv` (`readdlm(..., ',', skipstart = 1)`).
  - `const NDOTGV_NODAL_1D = "gx*pow(gx*gx + delta, -0.5)"`, `NDOTGV_NODAL_2D = "(g[1]*grad(v)[1] + g[2]*grad(v)[2])*pow(g[1]*g[1] + g[2]*g[2] + delta, -0.5)"`, `NDOTGV_ELEMENT_2D = "dot(grad(phi), grad(v))*pow(dot(grad(phi), grad(phi)) + delta, -0.5)"`.

- [ ] **Step 1: Write `tools/test_common.jl`** (uses `Test`; `include("../common.jl")` with `FINCH_TARGET` unset): `lumped_weights` on a 3×3 node grid with dx = 1 sums to 4.0 and the centre weight is 1.0; `finch_solver_kwargs(:lu)` returns `solver == :lu`; with `ENV["FINCH_TARGET"] = "mlir_cuda"` re-included in a fresh process (`run(\`julia +1.12 --project=$FINCH_V2_DIR -e ...\`)`) `finch_solver_kwargs(:lu).solver == :bicgstab`; `FINCH_TARGET=foo` in a fresh process exits non-zero and prints `FINCH_TARGET must be`.
- [ ] **Step 2: Run**: `julia +1.12 --project=$HOME/Packages/Finch-V2 tools/test_common.jl` — expected: fails, no `common.jl`.
- [ ] **Step 3: Implement `common.jl`** with the interfaces above.
- [ ] **Step 4: Run the test** — expected: all tests pass.
- [ ] **Step 5: Switch Example 1 to it**: in `main.jl` replace the target construction with `finch_target("poisson1d", PROJ_DIR)` and `init_finch("poisson1d"; finch_solver_kwargs(:lu)..., target = ...)`, `finch_run!(ctx, u)`; `run.sh` default `FINCH_V2_DIR="${FINCH_V2_DIR:-$HOME/Packages/Finch-V2}"` (already) and pass `FINCH_TARGET` through (it is inherited; add a usage comment `FINCH_TARGET=mlir_cuda ./run.sh`).
- [ ] **Step 6: Run and compare**: `cd Example1_Poisson1D && ./run.sh && python3 ../tools/compare_runlog.py ../regression/Example1_Poisson1D.txt output/runlog.txt` — expected: `OK`.
- [ ] **Step 7: Review-focus checks**: `FINCH_TARGET=foo ./run.sh` → non-zero, message lists the three values; `FINCH_TARGET=mlir_cpu ./run.sh` → non-zero, output contains `[build.sh]` (IREE missing on the Mac) and no Julia stack trace from post-processing.
- [ ] **Step 8: Commit**: `git add common.jl tools/test_common.jl Example1_Poisson1D && git commit -m "common.jl: FINCH_TARGET switch, run dispatch, shared post-processing; Example 1 uses it"`.

### Task 3: Examples 2, 3, 4 on the switch

**Files:**
- Modify: `Example2_Poisson2D/main.jl`, `Example3_ConvectionDiffusion1D/main.jl`, `Example4_ConvectionDiffusion2D/main.jl` (same three-line change as Example 1; Example 2 also renames `init_finch("heat2d")` to `"laplace2d"` and prints `max |u| = ...` and `mean u = ...` so it has numbers to compare from now on).

- [ ] **Step 1: Edit the three scripts** as in Task 2 step 5 (solver `:lu` kept).
- [ ] **Step 2: Run and compare** 3 and 4 against `regression/`: expected `OK`. Run 2; record its two printed numbers into `regression/Example2_Poisson2D.txt` (first reference for it).
- [ ] **Step 3: Commit**: `git add Example2_Poisson2D Example3_ConvectionDiffusion1D Example4_ConvectionDiffusion2D regression/Example2_Poisson2D.txt && git commit -m "Examples 2-4 on the FINCH_TARGET switch"`.

### Task 4: Example 5 (explicit CDI 1D) ported

**Files:**
- Modify: `Example5_ConservativeDiffusiveInterfaceExplicitTranslating1D/main.jl`

- [ ] **Step 1: Edit**: target switch (solver `:cg`); `time_stepper(STEP_EULER_EXPLICIT; Nsteps = NSTEPS, store_trajectory = true, snapshot_every = 25)` with `NSTEPS = round(Int, T / DT)` computed as before and `DT = T / NSTEPS` printed; snapshot/time-series writing through `read_snapshots` / `write_timeseries`; `NDOTGV_NODAL_1D` from common.
- [ ] **Step 2: Run and compare**: expected `OK` (`dt`, `max error`, `min phi`, `max phi`, snapshots).
- [ ] **Step 3: Commit**: `git commit -am "Example 5 on Finch V2 dev: Nsteps, common.jl"`.

### Task 5: Examples 6, 7, 8 ported one by one (no merge yet)

**Files:**
- Modify: `Example6_.../main.jl`, `Example7_.../main.jl`, `Example8_.../main.jl`

**Interfaces:**
- Produces the dev-form of the explicit 2D CDI script that Task 6 merges: unknowns `phi` and `g = variable("g"; kind = VAR_VECTOR, components = 2)`, `weak_form(g, "dot(g, v) - dot(grad(phi), v)")`, normal string `NDOTGV_NODAL_2D`, schedule `step_schedule(g, phi, monitor("mass", "phi"))` (element-wise branch: `step_schedule(phi, monitor(...))` with `NDOTGV_ELEMENT_2D`), `Nsteps`.

- [ ] **Step 1: Example 6**: apply the interface above (drop `gx, gy, ex, ey`), target switch, `Nsteps`, common post-processing. Run (≈10 min) and compare: expected `OK`.
- [ ] **Step 2: Commit** Example 6.
- [ ] **Step 3: Example 7**: same; keep `CDI_NORMAL`; `L1` via `lumped_weights`. Run (≈45 min, background, nothing else running) and compare: expected `OK`. Commit.
- [ ] **Step 4: Example 8**: same. Run (≈30 min) and compare: expected `OK`. Commit.
- [ ] **Step 5: If any comparison fails**: diff `output/finch/src/*.jl` of the new run against the old generated code (still in the old `output/finch` until the run overwrote it: copy `output/finch` to `output/finch_old` before step 1 of each example) to locate the changed operator, and report to the user before continuing.

### Task 6: Merge 6+7+8 into `Example6_ConservativeDiffusiveInterfaceExplicit2D`

**Files:**
- Rename: `Example6_ConservativeDiffusiveInterfaceExplicitCircle2D` → `Example6_ConservativeDiffusiveInterfaceExplicit2D` (`git mv`)
- Modify: its `main.jl` (case block), `run.sh` (usage comment), `matlabscripts/` (merge the three folders; file names prefixed by case where they clash; paths `output/<case>/`)
- Move: `Example7_.../paper512.sh`, `elementwise_normal.sh`, `Example8_.../paper512.sh` → `Example6_.../paper512_zalesak.sh`, `elementwise_normal.sh` (zalesak), `paper512_sheardrop.sh`, each exporting `CDI_CASE`
- Delete: `Example7_...`, `Example8_...` (`git rm -r`), after their regression references are renamed `regression/Example6_..._zalesak.txt`, `..._sheardrop.txt` (and the circle one `..._circle.txt`)

**Interfaces:**
- `CDI_CASE ∈ (circle, zalesak, sheardrop)`, default `circle`; unknown → `error("CDI_CASE must be circle, zalesak or sheardrop, got ...")`.
- Case block sets: `NEL` default (300 / 200 / 200), `T`, `IC(x, y)`, the velocity coefficient closures, `ADVECTION` string (circle: `"theta2*dot(u, grad(phi))*v"`; zalesak/sheardrop: `"- 0.5*theta2*phi*dot(u, grad(v)) - 0.5*theta2*phi*dot(u, grad(v))"`), `BCS` (circle: `boundary(phi, 1, BC_DIRICHLET, 0.0)`; others none), `snapshot_every` (40 / `max(1, NSTEPS ÷ 60)` / `max(1, NSTEPS ÷ 64)`), error reference (circle: translated tanh; zalesak/sheardrop: initial field, plus `phihalf` for sheardrop), `CDI_NORMAL` only honoured for zalesak.
- `OUT_DIR = joinpath(@__DIR__, get(ENV, "CDI_OUTDIR", joinpath("output", CASE)))`.

- [ ] **Step 1: `git mv` the folder**, write the case block, keep the printed summary lines byte-identical to the three old scripts (the comparer matches by line text).
- [ ] **Step 2: Run `CDI_CASE=circle`** and compare with `regression/Example6_..._circle.txt`: `OK`.
- [ ] **Step 3: Run `CDI_CASE=zalesak`** (45 min) and compare: `OK`. Check `output/circle/` still exists (review focus 3).
- [ ] **Step 4: Run `CDI_CASE=sheardrop`** (30 min) and compare: `OK`.
- [ ] **Step 5: `CDI_CASE=foo ./run.sh`** → non-zero with the case list.
- [ ] **Step 6: MATLAB paths**: `grep -n "output" Example6_*/matlabscripts/*.m` → every path is `../output/<case>/...`; no MATLAB run required.
- [ ] **Step 7: Commit**: `git add -A Example6_* regression && git rm -r Example7_* Example8_* && git commit -m "Example 6: explicit CDI 2D with cases circle, zalesak, sheardrop (merges old 6, 7, 8)"`.

### Task 7: Example 9 → Example 7 (semi-implicit CDI 1D) ported

**Files:**
- Rename: `Example9_ConservativeDiffusiveInterfaceImplicitTranslating1D` → `Example7_ConservativeDiffusiveInterfaceImplicit1D`
- Modify: `scheme.jl` (the scheme), `main.jl`, `convergence.jl`, `sharpening_sweep.jl` (only the include path and output names), `convergence.sh`, `sharpening_sweep.sh`, `run.sh`

**Interfaces:**
- `run_case(cfl, iters; form = "derived", nel = NEL, proj_dir, snapshots)` keeps its signature. Inside: fields `phi, phio, phis, gx` only; `phioo` → `old(phi, 2)` in the form; `a0, a1, a2` → coefficients `coefficient("a0", (x,y,z,t) -> ifelse(t < 1.5*dt_, 1/dt_, 1.5/dt_))` etc. with `dt_` the Julia constant (BDF1 on the first step, BDF2 after, exactly the old `update(a0 => ...)` sequence); `time_stepper(STEP_BDF2; Nsteps = nsteps, ...)` declared for the loop (forms contain no `Dt`, so they stay quasi-static; the stepper only supplies the loop and the history ring); `extrap` first lag `update(phis => "2*phio - old(phi, 2)")`; target switch with `:bicgstab`.

- [ ] **Step 1: `git mv`, edit `scheme.jl`** as above; `main.jl` unchanged except the include path.
- [ ] **Step 2: Run `./run.sh`** and compare with `regression/Example9_...txt` (renamed `regression/Example7_..._translating.txt`): `OK`. If `max error` differs in the 3rd digit, the `t` seen by the coefficient is not `t + dt` as assumed: switch the BDF switch to `ifelse(t < 0.5*dt_, ...)` and rerun; report which one reproduced.
- [ ] **Step 3: Smoke the drivers**: `julia +1.12 --project=$FINCH_V2_DIR sharpening_sweep.jl 100 derived 2 0.25` (one run, ≈10 s) writes `output/sharpening_sweep/results.csv` with one row; `convergence.jl` not rerun (hours), only `include`-checked by `julia -e 'include("convergence.jl")'` with `ARGS` that make it exit before running, if it has one; otherwise skip and say so in the commit message.
- [ ] **Step 4: Commit**: `git add -A Example7_* regression && git commit -m "Example 7: semi-implicit CDI 1D on dev (old Example 9): old(phi,2), BDF coefficients as t-dependent coefficients"`.

### Task 8: Examples 10, 11, 12 → Example 8 (semi-implicit CDI 2D), ported then merged

**Files:**
- Rename: `Example10_ConservativeDiffusiveInterfaceImplicitCircle2D` → `Example8_ConservativeDiffusiveInterfaceImplicit2D`
- Modify: `scheme.jl` (gains `case`), `main.jl`, `sharpening_sweep.jl`, `sharpening_sweep.sh`, `paper512.sh` (from old 12, `CDI_CASE=sheardrop`), `matlabscripts/` merged from 10, 11, 12
- Delete: `Example11_...`, `Example12_...`

**Interfaces:**
- `run_case(case::String, cfl, iters; form, nel, proj_dir, snapshots)`; `case ∈ (circle, zalesak, sheardrop)`; the case block as in Task 6 (advection string: circle `"theta2*dot(u, grad(phi))*v"`, others the split form; BCs; IC; velocity; T; snapshot cadence; error measure), scheme as in Task 7 with `g` vector variable (Task 5 interface) instead of `gx, gy, ex, ey`.
- `main.jl`: `CASE = get(ENV, "CDI_CASE", "circle")`, calls `run_case(CASE, 0.25, 2; ...)`, output `output/<case>/`.
- `sharpening_sweep.jl` passes `"circle"`.

- [ ] **Step 1: Port old 10 as is** (`git mv`, scheme edits, `g` vector), run `CDI_CASE=circle`, compare with the old 10 reference (renamed `regression/Example8_..._circle.txt`): `OK`. Commit.
- [ ] **Step 2: Add zalesak and sheardrop cases** from old 11 and 12 (their `main.jl` inline the same scheme; the case block carries the differences). Run each, compare with the renamed references: `OK`.
- [ ] **Step 3: `CDI_CASE=foo`** → non-zero with the case list. Smoke `sharpening_sweep.jl 100 derived 2 0.25`.
- [ ] **Step 4: MATLAB paths** `output/<case>/`.
- [ ] **Step 5: Commit**: `git add -A Example8_* regression && git rm -r Example11_* Example12_* && git commit -m "Example 8: semi-implicit CDI 2D with cases circle, zalesak, sheardrop (merges old 10, 11, 12)"`.

### Task 9: Examples 14, 15, 16 → 31 and 32 (CDI + dissolved gas), untracked

**Files:**
- Move (plain `mv`, folders are untracked): `Example14_...GasColumn1D` → `Example31_ConservativeDiffusiveInterfaceScalarTransportExplicit1D`; `Example15_...SeededBubbles2D` → `Example32_ConservativeDiffusiveInterfaceScalarTransportExplicit2D`; old 16 merged into 32 then deleted
- Create: `regression/Example16_...txt` from a baseline run (see step 1), then `regression/Example32_..._electrode.txt`
- Temporary: `/private/tmp/claude-501/.../scratchpad/Finch-V2-old` worktree of branch `mlir-batched-linalg` for the Example 16 baseline only, removed at the end of the task

**Interfaces:**
- Example 32: `CDI_CASE ∈ (bubbles, electrode)`, default `bubbles`. Case block: site list and IC, `CDI_AMIN` default (0.01 / 0.0), `CDI_T` default (1.5 / 2.0), mesh `periodic = :x` for electrode, BCs (all walls / `:ymax` only), the `psi` variable, its weak form and schedule entry only for electrode, advection terms `- alpha*dot(grad(psi), grad(v))` / `- c*dot(grad(psi), grad(v))` only for electrode, the electrode flux `- q*v*ds(:ymin)`, monitors (bubble regions / electrode input), `history.csv` columns, printed lines (kept byte-identical to old 15 and 16).
- Example 31 unchanged except `Nsteps`, switch, common helpers, `CDI_NEL` added (default 200).

- [ ] **Step 1: Example 16 baseline**: `git -C ~/Packages/Finch-V2 worktree add <scratchpad>/Finch-V2-old mlir-batched-linalg`; `FINCH_V2_DIR=<scratchpad>/Finch-V2-old ./run.sh` in old 16 (≈3 min); copy its summary lines to `regression/Example16_...txt`; `git -C ~/Packages/Finch-V2 worktree remove <scratchpad>/Finch-V2-old`.
- [ ] **Step 2: Port 14 → 31**, run, compare with the old 14 reference: `OK`.
- [ ] **Step 3: Port 15 → 32 (bubbles)**, run (≈2 min), compare: `OK`.
- [ ] **Step 4: Add the electrode case** from old 16, run (≈3 min), compare with the step-1 reference: `OK`. Run bubbles again: `OK` (the case block did not disturb it). `CDI_CASE=foo` → non-zero.
- [ ] **Step 5: Delete old 16's folder** (untracked: `rm -r`), MATLAB paths `output/<case>/`.
- [ ] **Step 6: No commit** (user rule). Note in the final report that 31/32 and their references are ready to commit when released.

### Task 10: Example 17 → Example 21 (Navier-Stokes + CDI)

**Files:**
- Rename: `Example17_NavierStokesConservativeDiffusiveInterfaceExplicit2D` → `Example21_NavierStokesConservativeDiffusiveInterfaceExplicit2D` (`git mv`; the untracked `output_case*` folders and `runNova.sh` move with it)
- Modify: `main.jl` (switch; `Nsteps`; common `read_snapshots`/`write_timeseries`; `dtc`/`idt` kept), `run.sh` (default `FINCH_V2_DIR=$HOME/Packages/Finch-V2`, cases A-E in the comment), `dump_snapshots.jl` (`read_snapshots`), `runNova.sh` (paths: `Example21_...`, `Finch-V2`; `JULIA_DEPOT_PATH` line kept)

- [ ] **Step 1: Edit**, run `CDI_CASE=A` (64², ≈3 min), compare with `regression/Example17_..._caseA.txt` (renamed `Example21_..._caseA.txt`): `OK`.
- [ ] **Step 2: Commit**: `git add -A Example21_* regression && git commit -m "Example 21: Navier-Stokes + CDI (old Example 17) on the FINCH_TARGET switch, Nsteps, common.jl"` (runNova.sh included; `output_case*` stay ignored).

### Task 11: Delete Example 13, README, clean.sh, final checks

**Files:**
- Delete: `Example13_InterfaceNormalFDvsFEM` (`git rm -r`)
- Modify: `README.md` (sections: Examples table with new numbers and cases; an "Old numbering" mapping table 1-17 → new; Setup with `git clone -b dev ... ~/Packages/Finch-V2`, commit 172f4fb, `Pkg.instantiate`; Run with `FINCH_TARGET` and `CDI_CASE`; Reference numbers pointing to `regression/`; the run-time paragraph updated; Nova pointer to `NOVA_SSH_SETUP_MAC.md` and "job scripts: see the next plan"), `clean.sh` (nothing to change if it works on ignored files; verify), `.gitignore` (add `output_old/`, `finch_old/` if used)
- Modify: `NovaGPU/runNovaGPU.sh` → delete the folder (the GPU test lives in Finch now; the next plan writes per-example GPU scripts)

- [ ] **Step 1: `git rm -r Example13_*`**; `rm -rf NovaGPU` (untracked locally; on Nova it is deleted in step 5).
- [ ] **Step 2: README** as above; the 31/32 rows are written but the README is committed only after the user releases 31/32 — so: write the README now, show the diff to the user, and ask in the final report whether the whole README may be committed.
- [ ] **Step 3: `./clean.sh -n`** lists only ignored files; `git status --short` shows only the expected untracked items (31, 32, their references, README).
- [ ] **Step 4: Commit** `git commit -m "Remove Example 13 (presentation-only); renumbering table"` with the deletion only; push after the user's answer on the README: `git push`.
- [ ] **Step 5: Nova**: `ssh nova-key 'cd /work/mech-ai/suresh/projects/FinchExamples && git pull && rm -rf NovaGPU && cd Example1_Poisson1D && FINCH_V2_DIR=/work/mech-ai/suresh/packages/Finch-V2 PATH=/work/mech-ai/suresh/packages/julia-1.12.7/bin:$PATH JULIA_DEPOT_PATH=/work/mech-ai/suresh/packages/julia-depot bash -c "julia --project=\$FINCH_V2_DIR main.jl"'` → prints `max error = 4.604168354349625e-7` (the Nova clone and depot work with 172f4fb; the Nova `run.sh` wrapper is the next plan's job, since `julia +1.12` is a juliaup call that Nova does not have).
