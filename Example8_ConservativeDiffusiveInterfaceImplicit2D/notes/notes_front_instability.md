# Example 8: the leading-edge wrinkle (diagnosis, 2026-10-16)

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
  Example 8 as derived                300^2   480   0.128       onset t~0.3
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
accuracy and the same second order in time (Examples 7 and 9, convergence.jl).
Not caused by: GMRES tolerance, linear solver, Picard count, time step, step count,
the outflow wall, the delta regularization of the normal, or tail wiggles.

## Setup pitfall
With R = 0.15 do not start the circle at xc0 <= 0.15: the back edge sits on the inflow
Dirichlet wall and 0.3% of the mass is cut off at the first step.
