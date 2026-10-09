# Finch, patched

This folder is a copy of [Finch](https://github.com/paralab/Finch) at upstream
commit `defa03c` (2023-06-06, "abs(detJ) in geometric factors") with a local
patch applied. It is kept inside this repository so that the examples run with
one clone and no separate Finch installation.

Patch date: 2026-09-25. Patch file: `fix-fem-boundary-integrals.patch`
(apply to a clean checkout of `defa03c` with `git apply`).

Finch is MIT licensed; see `LICENSE`. Nothing in this folder is original work
except the changes listed below.

## Why

A `boundary(...)` surface integral in an FEM `weakForm` never compiled on the
upstream code. Example 6 (2D CDI) uses one for its outflow boundary condition.

## What was changed (8 fixes, 3 files, +50 / -22 lines)

`src/IR_build_FEM.jl`

1. **:468** generated code read `geometric_factors.face_detj`; the field is `face_detJ`.
2. **:570, 574, 578, 582, 586, 590** `if !(x === nothing)` let an empty term
   vector through to `term_vec[1]` (BoundsError). Added `&& !is_empty_expression(x)`,
   as lines 394-401 and 475/482 already do for surface terms.
3. **:1027, 1030** `face2glb` was indexed by `:left_el` / `:right_el`; its third
   index is a face id everywhere else in Finch. Both changed to `:fid`.
4. **:1277, 1281, 1310, 1314, 1406** face loops emitted the element loop's bare
   `nodeID`; changed to `nodeID_L`, which the face loop defines.
5. after the `constantJ` block in `generate_term_calculation_fem` (~:1830): for
   `vors == "surface"`, use `face_wg` and `face_detj` instead of the volume
   `wg` and `detj`; and in the face loop body (~:470) add
   `face_wg = refel.surf_wg[face_refel_ind]` next to the existing `face_detj`.
6. **:1297-1307** the "no derivatives" surface-coefficient branch allocated its
   arrays with `nodes_per_face` while the loops fill `nodes_per_element` /
   `qnodes_per_element` entries; allocate at the loop sizes.

`src/generate_code_layer_julia.jl`

7. **:475** the `COEF_EVAL` fast path hardcoded `"(x,y,z,t,"`; it now uses
   `IR.args[4]`-`IR.args[7]`, as the `evaluate_coefficient` fallback below
   already did. (The GPU layer has the same hardcoding, untouched.)

`src/SymbolicParser.jl`

8. **~:466** "boundary integrals are not modified": boundary, Dirichlet-boundary
   and Neumann-boundary integral terms bypassed `reformat_for_stepper`, so they
   never received the `dt` factor (nor, for explicit steppers, the move of
   unknown-side terms to the right-hand side). They are now passed through in the
   surface slots with empty volume terms.

## Verification

`phi0 = x` on the unit square, `u = (1,0)`, pure advection with
`+ boundary(po*dot(u, normal())*v)`: one explicit step gives
`deltaM = -dt * (1 - 0.2113 h)` at every `dt` and every mesh (nel = 10..80, to
five digits). The exact answer is `-dt`.

## Known remaining issue (not fixed here)

The factor `(1 - 0.2113 h)` is structural: face integrands are evaluated at the
nearest *volume* Gauss point, via `face2local` into the volume basis matrix `Q`,
rather than on the face. Boundary integrals are therefore first-order accurate
regardless of element order, and a coefficient used *only* in a boundary term is
mis-evaluated (a constant 1 comes out as `(1 + 1/sqrt(3))/2 = 0.7887`), so a
spatially varying boundary mask cannot be used. The fix is to interpolate and
assemble face terms with `refel.surf_Q[face_refel_ind]` and place face-node values
at `face2local` slots; this is a restructuring of the face path and has been
reported to the Finch developers.

Also: a weak form whose only terms are boundary integrals (no volume terms) still
crashes. Degenerate case, not pursued.

## Manifest.toml

Upstream does not track `Manifest.toml`. It is included here so that
`Pkg.instantiate()` reproduces the dependency versions the examples were run with.
