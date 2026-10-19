# Example 13: run time (measured 2026-10-16)

300x300, CFL 0.25, 480 steps, two GMRES solves per step: ~8.4 s per step, about
65 minutes in total. Example 9 (same case without stabilization): 0.35 s per step.

Where the time goes (100x100 mesh, 40 steps, 80 solves):

    as written                          814 ms / solve
    tau replaced by a constant          777 ms / solve   -> tau's sqrt is not the cost
    direct solver instead of GMRES      811 ms / solve   -> the linear solve is not the cost

So the cost is the assembly of the stabilized weak form itself. Finch expands
each product symbolically, and tau * L(phi) * (u.grad v) with nhat written out
inline (it cannot be a parameter, see Example 12) produces a very large element
computation, roughly 20-70x the work per node of Example 9. Code generation also
grows to ~9 s.

Possible reductions, not done here because they change the formulation or its
fidelity to the derivation:
  - tau and/or nhat as nodal fields computed in Julia each Picard iteration and
    read through coefficient("...", "Main.array[node_index]"); a nodal (recovered)
    normal would then require the curvature term kappa* of the derivation.
  - merging the two stabilization halves, which coincide for a uniform velocity.

# Result (300x300, CFL 0.25, t = 0.4), measured 2026-10-16

    max error 0.0972     (Example 9 as derived: 0.128; Example 11 fully lagged: 0.0218;
                          Example 6 explicit: 0.036)
    phi in [0, 1] to 5e-10, relative mass change -3.1e-11

Leading edge, r(theta) - R in cells (the 3-element zigzag of Example 9):
    Ex 8  (no VMS):      0:+0.55  2:-0.16  4:+0.39  6:-0.21  8:+0.32  10:-0.20
    Ex 12 (VMS):         0:+0.39  2:-0.03  4:+0.30  6:-0.05  8:+0.19  10:+0.02
    Ex 10 (fully lagged) 0:+0.07  2:+0.07  4:+0.06  6:+0.05  8:+0.05  10:+0.05
VMS roughly halves the amplitude of the mode but does not remove it; the mode
belongs to the implicit factor of the sharpening term, which the stabilization
(built from the advection operator only) does not act on directly.
