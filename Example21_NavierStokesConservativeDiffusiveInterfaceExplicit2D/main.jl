#==============================================================================#
# 2D incompressible Navier-Stokes coupled to the Conservative Diffuse
# Interface (CDI) phase field, fully explicit, with an incremental projection
# method. Finch V2, development branch (budvinchathura/Finch-V2, branch dev),
# Julia target. Equations of the P30 note "Navier-Stokes Equations Coupled
# with the CDI Model, Explicit Coupling" (dimensional form, theta flags):
#
#   theta11 div(u) = 0
#   theta21 d(rho u)/dt + theta22 div(rho u u) - theta23 (rho1 - rho2) div(R u)
#     - theta24 div(mu grad u) - theta25 div(mu grad u^T) + theta26 grad p
#     - theta27 phi(1-phi)(1-2phi) grad phi + theta28 lap(phi) grad phi
#     - theta29 rho g = 0
#   theta31 dphi/dt + theta32 div(u phi) - div(R) = 0
#   R = theta33 gamma eps grad phi - theta34 gamma phi(1-phi) nhat
#
# rho = rho2 + (rho1 - rho2) phi, mu = mu2 + (mu1 - mu2) phi; phase 1 is
# phi = 1. The theta23 term moves momentum with the mass that the CDI flux
# R carries across the interface (Mirjalili and Mani, JCP 426 (2021) 109918),
# so the momentum flux is U u with U = rho u - (rho1 - rho2) R.
#
# One time step (the flowchart of the note), every entry a single linear solve:
#   0. alphao <- alpha                 keep phi^o (rho^o, mu^o, R^o, capillary)
#   1. g      projected gradient of phi^o (nodal normal, as in Examples 5-8)
#   2. alpha  CDI step with u^o (forward Euler, conservative advection)
#   3. ut     momentum predictor (Step 1): rho/dt ut = known terms, p~ = p^o
#   4. phat   pressure Poisson (Step 2): (dt/rho grad phat, grad q) = (ut, grad q)
#   5. u      velocity correction (Step 3): u = ut - dt/rho grad phat
#   6. p      <- p + phat
# The phase field is called alpha in the code (as in Examples 14-16).
#
# Cases (CDI_CASE):
#   A  drop advection (Mirjalili and Mani 2021, Case I, in 2D): a drop with
#      rho1/rho2 = 10 in a periodic box, uniform u = (1, 0), no viscosity or
#      surface tension. The exact solution keeps u = (1, 0) and the shape.
#      CDI_THETA23=0 gives the inconsistent momentum transport for comparison.
#   B  static drop: Laplace pressure jump sigma/R and spurious currents.
#   C  rising bubble, Hysing et al., IJNMF 60 (2009) 1259, test case 1:
#      [0,1]x[0,2], bubble R = 0.25 at (0.5, 0.5), rho 1000/100, mu 10/1,
#      g = 0.98, sigma = 24.5, no-slip top and bottom, free slip on the sides.
#   D  dam break (Martin and Moyce 1952, column of width a and height 2a):
#      box [0,4a]x[0,3a] with a = 0.25, water column in [0,a]x[0,2a], free
#      slip on all walls. Non-dimensional with the box width and g = 1:
#      rho 1/0.001, mu 1e-3/1.8e-6, sigma 2.2e-5 (water/air at 0.584 m, but
#      the water viscosity ten times the physical one, see the case block).
#      Front position and column height against the experiment, in
#      Z = x/a, H = h/(2a), T = t sqrt(2g/a).
#   E  falling drop: a water drop (R = 0.1) released at rest in air at
#      (0.5, 1.5) above a pool of depth 0.5 in [0,1]x[0,2], free slip walls,
#      the fluids of case D. Before impact the drop falls freely (v = g t,
#      air drag negligible), a check of the momentum coupling at density ratio
#      1000; after impact a qualitative splash. Bond number rho g D^2/sigma = 2.
#
# Pressure: all three cases have no pressure boundary condition (periodic or
# closed box), so the Poisson problem is singular (constants), and Finch V2
# does not handle that yet (CG stalls). A tiny reaction term preg*phat*q,
# 1e-6 of the Poisson diagonal (CDI_PREG), pins the constant.
#
# Continuous Galerkin, linear quadrilaterals, equal order for u and p.
# Run with ./run.sh (Julia 1.12, Finch V2 dev); FINCH_TARGET picks the CPU Julia or
# the MLIR CPU/GPU target (see ../common.jl). CDI_NEL and CDI_T override the
# mesh and the final time.
#==============================================================================#

include(joinpath(@__DIR__, "..", "common.jl"))   # Finch, FINCH_TARGET switch, helpers

const CASE = uppercase(get(ENV, "CDI_CASE", "A"))
CASE in ("A", "B", "C", "D", "E") || error("CDI_CASE must be A, B, C, D or E")

const OUT_DIR  = joinpath(@__DIR__, get(ENV, "CDI_OUTDIR", "output_case" * CASE))
const PROJ_DIR = joinpath(OUT_DIR, "finch")
mkpath(PROJ_DIR)

envf(k, d) = parse(Float64, get(ENV, k, string(d)))

# ---- case parameters ---------------------------------------------------------
# LY: domain height (width 1); NEL: elements across the width
if CASE == "A"
    const NEL = parse(Int, get(ENV, "CDI_NEL", "64"))
    const LY  = 1.0
    const RHO1, RHO2, MU1, MU2 = 10.0, 1.0, 0.0, 0.0
    const SIGMA, GRAV = 0.0, 0.0
    const XC, YC, RAD = 0.5, 0.5, 0.15        # D = 0.3 as in the paper
    const DROP_INSIDE = true                  # phase 1 (phi = 1) is the drop
    const U0  = (1.0, 0.0)
    const UREF = 1.0                          # velocity scale for gamma and dt
    const T_END = envf("CDI_T", 1.0)          # one period
elseif CASE == "B"
    const NEL = parse(Int, get(ENV, "CDI_NEL", "64"))
    const LY  = 1.0
    const RHO1, RHO2, MU1, MU2 = 10.0, 1.0, 0.01, 0.01
    const SIGMA, GRAV = 1.0, 0.0
    const XC, YC, RAD = 0.5, 0.5, 0.25
    const DROP_INSIDE = true
    const U0  = (0.0, 0.0)
    const UREF = 0.1
    const T_END = envf("CDI_T", 1.0)
elseif CASE == "D"
    const NEL = parse(Int, get(ENV, "CDI_NEL", "128"))
    const LY  = 0.75                           # 4a x 3a with a = 0.25
    const RHO1, RHO2 = 1.0, envf("CDI_RHO2", 1e-3)                  # water, air
    # water nu = 1e-3 (Re ~ 1500 on the box width): with the physical 1e-4 the
    # unstabilised Galerkin convection grows a mode at the air side of the face
    # (blows up at t ~ 0.2 for any dt and density ratio 10-1000); 1e-3 damps it
    const MU1, MU2   = envf("CDI_MU1", 1e-3), envf("CDI_MU2", 1.8e-6)
    const SIGMA, GRAV = 2.2e-5, 1.0
    const XC, YC, RAD = 0.25, 0.5, 0.0        # column width a = XC, height 2a = YC
    const DROP_INSIDE = true                  # phase 1 (phi = 1) is the water
    const U0  = (0.0, 0.0)
    const UREF = envf("CDI_UREF", 1.5)        # shallow-water front speed 2 sqrt(g 2a) = 1.4
    const T_END = envf("CDI_T", 2.0)
elseif CASE == "E"
    const NEL = parse(Int, get(ENV, "CDI_NEL", "96"))
    const LY  = 2.0
    const RHO1, RHO2 = 1.0, envf("CDI_RHO2", 1e-3)                  # water, air
    # water nu = 5e-3: with 1e-3 (case D) the fall is fine but the impact grows
    # the air-side mode and the run blows up at t ~ 1.46; 5e-3 carries the splash
    const MU1, MU2   = envf("CDI_MU1", 5e-3), envf("CDI_MU2", 1.8e-6)
    const SIGMA, GRAV = 0.02, 1.0             # Bo = rho g (2R)^2 / sigma = 2
    const XC, YC, RAD = 0.5, 1.5, 0.1
    const POOL = 0.5                          # pool surface
    const DROP_INSIDE = true
    const U0  = (0.0, 0.0)
    const UREF = envf("CDI_UREF", 1.5)        # free fall over 0.9: sqrt(2 g 0.9) = 1.34
    const T_END = envf("CDI_T", 2.5)
else
    const NEL = parse(Int, get(ENV, "CDI_NEL", "80"))
    const LY  = 2.0
    const RHO1, RHO2, MU1, MU2 = 1000.0, 100.0, 10.0, 1.0   # liquid, gas
    const SIGMA, GRAV = 24.5, 0.98
    const XC, YC, RAD = 0.5, 0.5, 0.25
    const DROP_INSIDE = false                 # phase 1 = liquid outside the bubble
    const U0  = (0.0, 0.0)
    const UREF = 0.3                          # rise velocity peaks at 0.24
    const T_END = envf("CDI_T", 3.0)
end

const DX   = 1.0 / NEL
const EPS  = envf("CDI_EPS_DX", 1.0) * DX     # interface thickness
const GAM  = envf("CDI_GAMMA", 1.0) * UREF    # sharpening velocity
const TH23 = envf("CDI_THETA23", 1.0)         # 0: inconsistent momentum flux
# Pressure pin: preg*phat*q with preg = PREG_REL times the Poisson diagonal
# scale dt/(rho_max dx^2), so it fixes the constant and nothing else
const PREG_REL = envf("CDI_PREG", 1e-6)

# Time step: the smallest of the four limits of the note, with 0.2 for the
# consistent mass matrix on the phase field limit and 0.25 on the others
lims = Dict("phase field" => 0.2 * DX^2 / (4 * GAM * EPS),
            "convective"  => 0.25 * DX / UREF)
# Viscous: forward Euler on M^-1 K, consistent mass (max eigenvalue 24/dx^2
# for Q1 in 2D, three times the lumped one) and the full stress mu(grad u +
# grad u^T) (twice the Laplacian's on compressive modes): dt <= dx^2/(24 nu).
# The note's dx^2 rho/mu with C = 1/16 was unstable at t ~ 0.07 in case C.
MU1 + MU2 > 0  && (lims["viscous"]   = 0.7 * DX^2 * minimum(r / m for (r, m) in ((RHO1, MU1), (RHO2, MU2)) if m > 0) / 24)
SIGMA > 0      && (lims["capillary"] = 0.25 * sqrt((RHO1 + RHO2) * DX^3 / (2pi * SIGMA)))
const DT_LIM = envf("CDI_DTFAC", 1.0) * minimum(values(lims))   # CDI_DTFAC scales dt
const NSTEPS = ceil(Int, T_END / DT_LIM - 1e-9)
const DT     = T_END / NSTEPS

# ---- Finch set-up --------------------------------------------------------------
# Every matrix is SPD (density-weighted mass, mass, variable Poisson): CG
# CDI_MAXIT: CG iteration cap (the Jacobi-preconditioned pressure solve needs
# more than Finch's default 500 at 120 x 240)
ctx = init_finch("nscdi2d"; verbose = false,
                 finch_solver_kwargs(:cg; tol = 1e-10, max_iter = parse(Int, get(ENV, "CDI_MAXIT", "5000")))...,
                 target = finch_target("nscdi", PROJ_DIR))

domain(2)
if CASE in ("C", "D", "E")
    m = mesh(MESH_QUAD; elsperdim = (NEL, round(Int, LY * NEL)))
    m.nodes[2, :] .*= LY       # the built-in mesh is the unit square: stretch it
else
    m = mesh(MESH_QUAD; elsperdim = NEL, periodic = (:x, :y))
end
function_space(order = 1)

alpha  = variable("alpha")                     # phase field phi
alphao = variable("alphao")                    # phi^o, set by update
g      = variable("g"; kind = VAR_VECTOR)      # projected grad phi^o
u      = variable("u"; kind = VAR_VECTOR)      # velocity
ut     = variable("ut"; kind = VAR_VECTOR)     # predicted velocity
phat   = variable("phat")                      # pressure correction
p      = variable("p")                         # pressure (auxiliary, update)
test_symbol("v")

if CASE == "C"
    for b in (:ymin, :ymax)                    # no slip
        boundary(u,  b, BC_DIRICHLET, 0.0)
        boundary(ut, b, BC_DIRICHLET, 0.0)
    end
    for b in (:xmin, :xmax)                    # free slip: u_x = 0, no shear
        boundary(u[1],  b, BC_DIRICHLET, 0.0)
        boundary(ut[1], b, BC_DIRICHLET, 0.0)
    end
elseif CASE in ("D", "E")                      # free slip on all four walls
    for b in (:xmin, :xmax)
        boundary(u[1],  b, BC_DIRICHLET, 0.0)
        boundary(ut[1], b, BC_DIRICHLET, 0.0)
    end
    for b in (:ymin, :ymax)
        boundary(u[2],  b, BC_DIRICHLET, 0.0)
        boundary(ut[2], b, BC_DIRICHLET, 0.0)
    end
end

sgn = DROP_INSIDE ? 1.0 : -1.0
if CASE == "D"
    # water column [0, XC] x [0, YC]: signed distance to its boundary (exact
    # inside and beside the faces, slightly short outside the corner)
    phi0(x, y) = 0.5 * (1 - tanh(max(x - XC, y - YC) / (2 * EPS)))
elseif CASE == "E"
    # drop and pool, both phase 1: the larger of the two profiles
    phi0(x, y) = max(0.5 * (1 - tanh((sqrt((x - XC)^2 + (y - YC)^2) - RAD) / (2 * EPS))),
                     0.5 * (1 - tanh((y - POOL) / (2 * EPS))))
else
    phi0(x, y) = 0.5 * (1 - sgn * tanh((sqrt((x - XC)^2 + (y - YC)^2) - RAD) / (2 * EPS)))
end
initial(alpha,  (x, y, z, t) -> phi0(x, y))
initial(alphao, (x, y, z, t) -> phi0(x, y))
initial(u, [(x, y, z, t) -> U0[1], (x, y, z, t) -> U0[2]])

# theta flags, dimensional form
coefficient("th11", 1.0); coefficient("th21", 1.0); coefficient("th22", 1.0)
coefficient("th23", TH23); coefficient("th24", 1.0); coefficient("th25", 1.0)
coefficient("th26", 1.0); coefficient("th29", 1.0)
coefficient("th27", SIGMA > 0 ? 6 * SIGMA / EPS : 0.0)
coefficient("th28", 6 * SIGMA * EPS)
coefficient("th31", 1.0); coefficient("th32", 1.0)
coefficient("th33", 1.0); coefficient("th34", 1.0)
coefficient("rho1", RHO1); coefficient("rho2", RHO2)
coefficient("mu1", MU1);   coefficient("mu2", MU2)
coefficient("gam", GAM);   coefficient("eps", EPS)
coefficient("delta", 1e-12)
coefficient("preg", PREG_REL * DT / (max(RHO1, RHO2) * DX^2))
# dt is reserved in the form language but is not the step in a quasi-static
# solve, so the step size enters as coefficients
coefficient("dtc", DT); coefficient("idt", 1 / DT)
coefficient("grav", [(x, y, z, t) -> 0.0, (x, y, z, t) -> -GRAV])

time_interval(T_END)
const NFRAMES = parse(Int, get(ENV, "CDI_FRAMES", "60"))
time_stepper(STEP_EULER_EXPLICIT; Nsteps = NSTEPS, store_trajectory = true,
             snapshot_every = max(1, NSTEPS ÷ NFRAMES))

# ---- weak forms ------------------------------------------------------------------
# Index convention: grad(w)[i, j] = d w_i / d x_j
IJ = [(i, j) for i in 1:2 for j in 1:2]
sumij(f) = "(" * join((f(i, j) for (i, j) in IJ), " + ") * ")"

const RHO  = "(rho2 + (rho1 - rho2)*alpha)"     # current level (after step 2)
const RHOO = "(rho2 + (rho1 - rho2)*alphao)"
const MUO  = "(mu2 + (mu1 - mu2)*alphao)"
const INVG = "pow(g[1]*g[1] + g[2]*g[2] + delta, -0.5)"
nhat(j) = "g[$j]*" * INVG
# CDI flux R^o, with the same normal the phase field step uses
Ro(j) = "(th33*gam*eps*grad(alphao)[$j] - th34*gam*alphao*(1 - alphao)*" * nhat(j) * ")"
const DIVV = "(grad(v)[1, 1] + grad(v)[2, 2])"
# grad(phi^o) in the capillary stress: the L2-projected nodal gradient g
# (proj, default) or the element gradient (elem). With elem the static drop
# (case B) blows up by t ~ 0.05 whatever dt; with proj it is stable. The
# diffusion part of R stays element-wise in both: it must be the compact
# Laplacian of the phase field step, which damps grid-scale modes, and R^o
# must equal that step's flux.
const STGRAD = get(ENV, "CDI_STGRAD", "proj")
STGRAD in ("elem", "proj") || error("CDI_STGRAD must be elem or proj")
gphi(i) = STGRAD == "proj" ? "g[$i]" : "grad(alphao)[$i]"

# 1. nodal gradient of phi^o
weak_form(g, "dot(g, v) - dot(grad(alphao), v)")

# 2. phase field (alphao = alpha here, at the previous level)
weak_form(alpha, "th31*Dt(alpha)*v - th32*alpha*dot(u, grad(v))" *
                 " + th33*gam*eps*dot(grad(alpha), grad(v))" *
                 " - th34*gam*(alpha - alpha*alpha)*(g[1]*grad(v)[1] + g[2]*grad(v)[2])*" * INVG)

# 3. momentum predictor, Eq. (step1) of the note moved to one side
weak_form(ut, "th21*$RHO*idt*dot(ut, v) - th21*$RHOO*idt*dot(u, v)" *
    " - th22*$RHOO*" * sumij((i, j) -> "u[$j]*u[$i]*grad(v)[$i, $j]") *
    " + th23*(rho1 - rho2)*" * sumij((i, j) -> Ro(j) * "*u[$i]*grad(v)[$i, $j]") *
    " + th24*$MUO*" * sumij((i, j) -> "grad(u)[$i, $j]*grad(v)[$i, $j]") *
    " + th25*$MUO*" * sumij((i, j) -> "grad(u)[$j, $i]*grad(v)[$i, $j]") *
    " - th26*p*$DIVV" *
    " + th27*0.5*alphao*alphao*(1 - alphao)*(1 - alphao)*$DIVV" *
    " - th28*" * sumij((i, j) -> gphi(j) * "*" * gphi(i) * "*grad(v)[$i, $j]") *
    " + th28*0.5*(" * gphi(1) * "*" * gphi(1) * " + " * gphi(2) * "*" * gphi(2) * ")*$DIVV" *
    " - th29*$RHOO*dot(grav, v)")

# 4. pressure Poisson for the correction; walls and periodic sides: g^p = 0
weak_form(phat, "th11*th26*dtc/(th21*$RHO)*dot(grad(phat), grad(v))" *
                " - th11*dot(ut, grad(v)) + preg*phat*v")

# 5. velocity correction, an L2 projection
weak_form(u, "dot(u, v) - dot(ut, v) + th26*dtc/(th21*$RHO)*dot(grad(phat), v)")

# Bubble (case C) or drop (A, B) indicator for the monitors
const IND = DROP_INSIDE ? "alpha" : "(1 - alpha)"

step_schedule(update(alphao => "alpha"), g, alpha, ut, phat, u,
              update(p => "p + phat"),
              monitor("volume",   "alpha"),
              monitor("kinetic",  "0.5*$RHO*dot(u, u)"),
              monitor("div2",     "(grad(u)[1, 1] + grad(u)[2, 2])*(grad(u)[1, 1] + grad(u)[2, 2])"),
              monitor("ind",      IND),
              monitor("ind_x",    "x*$IND"),
              monitor("ind_y",    "y*$IND"),
              monitor("ind_ux",   "u[1]*$IND"),
              monitor("ind_uy",   "u[2]*$IND"),
              monitor("perim",    "pow(dot(grad(alpha), grad(alpha)) + delta, 0.5)"),
              monitor("uerr2",    "(u[1] - $(U0[1]))*(u[1] - $(U0[1])) + (u[2] - $(U0[2]))*(u[2] - $(U0[2]))"))

println("case $CASE: NEL = $NEL, LY = $LY, dt = $DT, steps = $NSTEPS, theta23 = $TH23, preg_rel = $PREG_REL, stgrad = $STGRAD")
println("  dt limits: ", join(("$k $(round(v, sigdigits = 3))" for (k, v) in lims), ", "))
t0 = time()
solve()
finch_run!(ctx, alpha, u, p)
println("run time = ", round(time() - t0, digits = 1), " s")

# ---- output ---------------------------------------------------------------------
x = ctx.grid.nodes[1, :]; y = ctx.grid.nodes[2, :]; N = length(x)
results = joinpath(PROJ_DIR, "results")
times, fields = read_snapshots(PROJ_DIR)
check_final_time(times, T_END)
A  = fields["alpha"]
UX = fields["u"][1:N, :]          # vector fields are component-blocked
UY = fields["u"][N+1:2N, :]
P  = fields["p"]
write_timeseries(OUT_DIR, "alpha", [x y], A)
write_timeseries(OUT_DIR, "ux", [x y], UX)
write_timeseries(OUT_DIR, "uy", [x y], UY)
write_timeseries(OUT_DIR, "p", [x y], P)
write_snapshot_times(OUT_DIR, times)

# history.csv: t, volume, kinetic energy, ||div u||_2, centroid x, y,
# velocity x, y, circularity, ||u - u0||_2
mon = readdlm(joinpath(results, "monitors.csv"), ',', skipstart = 1)
t_m = mon[:, 2]
col(k) = mon[:, 2 + k]
vol, ke, div2, ind, ix, iy, iux, iuy, per, uerr2 = (col(k) for k in 1:10)
circ = 2 .* sqrt.(pi .* ind) ./ per
hist = [t_m vol ke sqrt.(div2) ix ./ ind iy ./ ind iux ./ ind iuy ./ ind circ sqrt.(uerr2)]
writedlm(joinpath(OUT_DIR, "history.csv"), hist, ',')

a = A[:, end]
println("volume: ", vol[1], " -> ", vol[end], ", relative change ", (vol[end] - vol[1]) / vol[1])
println("kinetic energy: ", ke[1], " -> ", ke[end], ", relative change ", (ke[end] - ke[1]) / max(ke[1], 1e-300))
println("min alpha = ", minimum(a), ", max alpha = ", maximum(a))
println("max ||div u||_2 = ", maximum(sqrt.(div2)))
println("max |u| at the end = ", maximum(sqrt.(UX[:, end] .^ 2 .+ UY[:, end] .^ 2)))
if CASE == "A"
    println("||u - (1,0)||_2 at the end = ", sqrt(uerr2[end]),
            ", max |u - (1,0)| = ", maximum(sqrt.((UX[:, end] .- 1) .^ 2 .+ UY[:, end] .^ 2)))
    l1 = sum(abs.(a .- A[:, 1])) * DX^2
    println("L1 shape error after one period = ", l1)
elseif CASE == "B"
    pin  = sum(P[k, end] for k in 1:N if hypot(x[k] - XC, y[k] - YC) < 0.5RAD) / count(k -> hypot(x[k] - XC, y[k] - YC) < 0.5RAD, 1:N)
    pout = sum(P[k, end] for k in 1:N if hypot(x[k] - XC, y[k] - YC) > 1.6RAD) / count(k -> hypot(x[k] - XC, y[k] - YC) > 1.6RAD, 1:N)
    println("pressure jump = ", pin - pout, " (Laplace sigma/R = ", SIGMA / RAD, ")")
elseif CASE == "C"
    k = argmax(hist[:, 8])
    kc = argmin(circ)
    println("max rise velocity = ", hist[k, 8], " at t = ", t_m[k], "  (Hysing TC1: 0.2417 at 0.9213)")
    println("min circularity   = ", circ[kc], " at t = ", t_m[kc], "  (Hysing TC1: 0.9013 at 1.9041)")
    println("centroid y at end = ", hist[end, 6], "  (Hysing TC1 at t = 3: 1.0813)")
end
if CASE == "E"
    # drop.csv: t, centroid y and mean vertical velocity of the water above the
    # pool (y > POOL + 4 eps), valid until the drop touches the pool; free fall
    # gives y_c = YC - g t^2/2, v = -g t
    above = findall(k -> y[k] > POOL + 4EPS, 1:N)
    yc = [sum(A[above, j] .* y[above]) / max(sum(A[above, j]), 1e-300) for j in 1:length(times)]
    vc = [sum(A[above, j] .* UY[above, j]) / max(sum(A[above, j]), 1e-300) for j in 1:length(times)]
    writedlm(joinpath(OUT_DIR, "drop.csv"), [times yc vc], ',')
    k = findfirst(j -> minimum(y[k] for k in 1:N if A[k, j] > 0.5 && y[k] > POOL + 4EPS; init = Inf) < POOL + 6EPS, 1:length(times))
    println("drop reaches the pool at t = ", k === nothing ? "never" : times[k], " (free fall: ", sqrt(2 * (YC - RAD - POOL) / GRAV), ")")
end
if CASE == "D"
    # front.csv: t, x of the front on the floor, water height on the left wall
    # (phi = 0.5 crossings, linearly interpolated), and the same in the
    # Martin-Moyce variables T = t sqrt(2g/a), Z = x/a, H = h/(2a)
    a = XC
    function crossing(s, v)          # last s where v crosses 0.5 downward along s
        o = sortperm(s); s = s[o]; v = v[o]; c = v[end] >= 0.5 ? s[end] : s[1]   # wet to the end: the wall
        for k in 1:length(s)-1
            if v[k] >= 0.5 && v[k+1] < 0.5
                c = s[k] + (v[k] - 0.5) / (v[k] - v[k+1]) * (s[k+1] - s[k])
            end
        end
        c
    end
    bottom = findall(k -> y[k] < 0.5 * DX, 1:N); left = findall(k -> x[k] < 0.5 * DX, 1:N)
    xf = [crossing(x[bottom], A[bottom, j]) for j in 1:length(times)]
    hl = [crossing(y[left], A[left, j]) for j in 1:length(times)]
    writedlm(joinpath(OUT_DIR, "front.csv"), [times xf hl times .* sqrt(2 * GRAV / a) xf ./ a hl ./ (2a)], ',')
    k = findfirst(z -> z >= 4a - 2DX, xf)
    println("front reaches the wall at t = ", k === nothing ? "never" : times[k], " (T = ", k === nothing ? "-" : times[k] * sqrt(2 * GRAV / a), "; experiment T ~ 3)")
end
println(length(times), " snapshots written")
