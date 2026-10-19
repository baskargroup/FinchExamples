# Example 9: the leading-edge wrinkle (diagnosis, 2026-10-16)

Error = max nodal |phi - phi_exact| at t = 0.4 against the translated tanh circle.

## What it is
A sawtooth along the interface, wavelength 3 elements, confined to the arc where the
front is parallel to the grid (|theta| < ~12 deg at the leading edge). Amplitude
doubles every ~0.04 time units from round-off: an instability of the scheme, not a
truncation error. r(theta) - R at t = 0.4 in cells:
  0:+0.55  2:-0.16  4:+0.39  6:-0.21  8:+0.32  10:-0.20  15:~0  180:+0.09 (smooth)
Centroid exact to 1e-6, shape exact to 1e-4, phi in [0,1], mass to 1e-11. Invisible
in whole-domain plots; see output/Figure5_FrontZoom.pdf.

## Runs (CFL 0.25, two Picard iterations unless stated)
  run                                 mesh   steps  max error   wrinkle
  explicit (Example 6)                300^2  2400   0.036       none
  Example 9 as derived                300^2   480   0.128       onset t~0.3
  3 Picard iterations                 300^2   480   0.125       unchanged
  CFL 0.125                           300^2   960   0.138       same onset TIME
  circle started at x = 0.2           300^2   480   0.128       identical (not the wall)
  delta 1e-8 / delta 1                300^2   480   0.088/0.111 smaller, still there
  BDF1 throughout                     300^2   480   0.594       far worse, onset t~0.13
  150^2 GMRES / direct                150^2   240   0.021       none (solver irrelevant)
  150^2, CFL 0.125                    150^2   480   0.021       none to t = 0.4
  150^2 to t = 0.6 (direct)           150^2   360   -           appears at t~0.45
  200^2                               200^2   320   0.0216      none to t = 0.38
  400^2                               400^2   640   0.293       onset t~0.2
  sharpening fully lagged (Ex. 10)    300^2   480   0.0218      NONE

## Cause
The semi-implicit split of the sharpening term, gamma*phi*(1-phis)*nhat(phis) with phi
implicit, coupled to the lagged normal, is anti-diffusive for 3-cell modes tangential
to a grid-aligned front. Evaluating the whole term from the Picard iterate,
gamma*phis*(1-phis)*nhat(phis), removes the mode; in 1D the two forms give the same
accuracy and the same second order in time (Examples 8 and 10, convergence.jl).
Not caused by: GMRES tolerance, linear solver, Picard count, time step, step count,
the outflow wall, the delta regularization of the normal, or tail wiggles.

## Setup pitfall
With R = 0.15 do not start the circle at xc0 <= 0.15: the back edge sits on the inflow
Dirichlet wall and 0.3% of the mass is cut off at the first step.

## Outflow bookkeeping (2026-10-16, open question)
Integrating Example 9's wall values, int_0^0.4 int phi(1,y) dy dt = 5.3e-13 absolute =
7.5e-12 relative (mass0 = 0.0708). Measured relative mass change at t = 0.4: +5.4e-12,
so solver drift (GMRES reltol 1e-12, 960 solves) ~ +1.3e-11 and outflow ~ -7.5e-12 are
the same order; the mass plot resolves neither. Example 6 (explicit, same wall values
~1.4e-10 at t = 0.4) measured -8.4e-13 relative for the same outflow term - a factor ~9
below this integral. Unexplained; candidates: the first-order face quadrature (bug 9,
face values interpolated through the volume basis), or the lagged g. Worth a direct
test (one step, known wall profile) before quoting either number as "the outflow".

## Iteration-count check across linearizations (2026-10-16)

The fully lagged form (Example 11) run with four Picard iterations instead of two gives
max error 0.067 with a zigzag at the TRAILING edge (180 deg: +0.28 el, neighbours
-0.06); the Newton form (Example 15) gives 0.055 with two iterations and 0.113 with
three, also at the trailing edge; the derived form gives 0.125 with three (leading
edge). So the grid-aligned mode is a property of the iterated semi-implicit step with
the normal rebuilt from the iterate; the fully lagged two-iteration result (0.022)
stops before it grows. See Example15/notes/notes_newton_2d.md.

## The explicit scheme has the same mode (2026-10-16, scratchpad/ex6rk)

Example 6 (explicit, 300 x 300, same case) with RK4 instead of forward Euler at the
same dt: max error 0.112 (Euler 0.036), zigzag at the leading edge (0:-0.39 4:+0.28
6:+0.29) AND the trailing edge (176:+0.28 180:-0.34 184:+0.28); mass -7e-12, bounds ok.
In 1D (Example 5) RK4 is simply better: 0.0047 vs 0.0156, stable to 2x the Euler dt.
Conclusion: the grid-aligned wave is a property of the Galerkin space discretization
of the CDI equation (sharpening with element-wise normal); forward Euler's first-order
error damps it, BDF2 with two under-converged Picard iterations of the fully lagged
form stops before it grows, and every accurate time integration (RK4, converged
Picard, Newton+delta with 4 iterations) shows it. It is not a time-scheme problem.

## FIX FOUND: recovered nodal normal (2026-10-17, scratchpad/ex8_nodalnormal)

Same scheme as Example 9 (BDF2, 2 Picard, derived sharpening split, GMRES 1e-12,
300 x 300, CFL 0.25, delta 1e-12) but nhat is a VECTOR variable computed at the nodes
by central differences of phis on the uniform grid (one-sided at the walls), refreshed
before each of the two solves, instead of the element-wise P1 gradient:

    max error 0.00675   (element-wise normal: 0.128; fully lagged 0.0218; explicit 0.036)
    front r - R in elements: leading 0..20 deg all within +-0.01, trailing 170..190 deg
    +0.01..+0.02, 45/90/135 deg +-0.01. No zigzag anywhere.
    phi in [0, 1 - 6e-10], mass -4.6e-11, run time 235 s.

So the grid-aligned mode comes from the element-wise normal (constant per element,
discontinuous across edges, noise in the tails), not from the time scheme or the
linearization. Mirjalili et al. (JCP 2020) compute the normal at cell centres by central
differences, i.e. a continuous field, which is what this reproduces. Next: same change
in the explicit Zalesak case (scratchpad/ex15_nodalnormal).

## Nodal normal in the explicit examples (2026-10-17; numbering is the NEW one)

Examples 5, 6 and 7 now recover the normal at the nodes by a lumped L2 projection of
the FE gradient (recover_nhat! in each main.jl, built from Finch's loc2glb, refel and
geo_factors; the finite-difference version used for the first tests gave the same
numbers to 2%).
  Example 5 (1D):     max error 0.015621 with either normal (identical to 9 digits:
                      in 1D the normal is +-1 wherever the gradient is nonzero).
  Example 6 (circle): max error 0.0183 (element-wise 0.0360); front r - R = -0.07 el
                      at every angle, i.e. a uniform small offset and no zigzag;
                      mass -5e-12 (CG default tolerance), bounds ok.
  Example 7 (Zalesak): L1 3.39e-3 (element-wise 1.18e-2), see its notes.
The semi-implicit Examples 8-16 still use the element-wise normal (pending).

## Examples 8 and 9 converted (2026-10-17): recover_nhat! from phis before each solve

  Example 8 (1D): max error 0.0058806 (unchanged to 7 digits); temporal order 2.32 with
                  two Picard iterations (2.27 before; only the CFL 1.0 point moved,
                  0.0465 -> 0.0563), 0.98 with one.
  Example 9 (2D): max error 0.00622 (element-wise normal 0.128; finite-difference nodal
                  normal 0.00675), mass -4.6e-11, bounds ok, 167 s. Front flat to
                  +-0.01 element at every angle.
Folders renamed ...ImplicitTranslating1D and ...ImplicitCircle2D. Examples 10-16 still
element-wise.
