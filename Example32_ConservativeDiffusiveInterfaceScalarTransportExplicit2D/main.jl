#==============================================================================#
# 2D CDI phase field with dissolved-gas transport, solved explicitly with
# Finch V2 (branch dev). The formulation is the P30 note "CDI with Scalar
# Transport, Explicit Galerkin FE Formulation for Diffusive Bubble Growth":
#
#   alpha_t + div(u*alpha) - div(gamma*eps*grad alpha)
#     + div(gamma*(alpha - alpha^2)*nhat) - theta15*J = 0
#
#   c_t + div(u*c) - div(D*grad c) - div(D*alpha*c/eps*nhat) + J = 0
#
#   J = alpha/tau * max(0, c - csat*(1 - alpha)),     div(u) = theta15*J
#
# Two cases (CDI_CASE):
#
#   bubbles   (default) growth of five seeded bubbles in a liquid supersaturated
#             at c = 3*csat, alpha = 0 and c = csuper on the whole boundary.
#             u = 0 (the note's choice): the bubbles grow in place by the volume
#             source alone, no liquid is displaced and the dissolved gas is
#             replenished only by diffusion. Three bubbles are isolated (radii
#             0.08, 0.06, 0.04); two of radius 0.06 sit close together and
#             compete for the gas between them. D is small so that the depletion
#             layers (~ sqrt(D*T) = 0.1) stay clear of the boundary; tau = 0.5.
#             T = 1.5 is about as far as this case goes: with u = 0 the dissolved
#             gas piles up around the bubbles, alpha exceeds 1 from t ~ 1.6 and
#             the explicit step blows up at t ~ 2.2 (100^2).
#
#   electrode bubbles growing on a gas-evolving electrode (after Qiu and Mani):
#             the bottom wall produces dissolved gas at a constant rate q per
#             unit area into a liquid that starts saturated (c = csat); the
#             bubbles seeded on the wall (half discs of different sizes) take
#             up the gas once the liquid next to the wall is supersaturated,
#             and close neighbours merge. Velocity: the expansion flow
#             u = grad(psi), lap(psi) = theta15*J, psi = 0 at the top and
#             d(psi)/dn = 0 on the wall (the flow around a growing bubble
#             without buoyancy or surface tension). Walls: bottom electrode
#             (flux q on c, zero flux on alpha, no flow through), top saturated
#             liquid (alpha = 0, c = csat, the displaced liquid leaves here),
#             sides periodic (a row of sites on a long electrode).
#
# One forward Euler step is a step_schedule, every entry reading the previous level:
#   1. update J    node by node, with the nucleation threshold alpha > AMIN
#   2. solve  psi  (electrode only, quasi-static) for the velocity
#   3. solve  g    the nodal gradient of alpha by L2 projection (consistent
#                  mass); nhat = g/sqrt(g.g + delta). The element-wise normal
#                  gives a grid-scale zigzag on interfaces along the grid.
#   4. step   c
#   5. step   alpha
#
# Trace gas in a supersaturated liquid grows like exp(2*theta15*t/tau) since J
# is proportional to alpha; the threshold AMIN (CDI_AMIN) removes that source.
# "max alpha away from the bubbles" / "in the upper liquid" is printed as a check.
#
# Continuous Galerkin, linear quadrilaterals. Run with ./run.sh (Julia 1.12);
# the target (CPU Julia, or MLIR on the CPU/GPU) is chosen by FINCH_TARGET, see
# ../common.jl. CDI_NEL and CDI_T override the mesh and the final time, CDI_AMIN
# the threshold, CDI_OUTDIR the output folder (default output/<case>/).
#==============================================================================#

include(joinpath(@__DIR__, "..", "common.jl"))   # Finch, FINCH_TARGET switch, helpers

const CASE = get(ENV, "CDI_CASE", "bubbles")
CASE in ("bubbles", "electrode") || error("CDI_CASE must be bubbles or electrode, got \"$CASE\"")
const ELECTRODE = CASE == "electrode"

# CSV results go to output/<case>/; Finch's generated project (and its raw
# results) to its finch/ subfolder
const OUT_DIR  = joinpath(@__DIR__, get(ENV, "CDI_OUTDIR", joinpath("output", CASE)))
const PROJ_DIR = joinpath(OUT_DIR, "finch")
mkpath(PROJ_DIR)

# Parameters (non-dimensional: domain size 1, concentration scaled by csat)
const NEL    = parse(Int, get(ENV, "CDI_NEL", "200"))
const L      = 1.0
const DX     = L / NEL
const CSAT   = 1.0
const CSUPER = 3.0 * CSAT       # bubbles: the supersaturated liquid
const Q      = 0.1              # electrode: dissolved gas produced, per unit area
const D      = ELECTRODE ? 0.005 : 0.01   # diffusivity of dissolved gas in the liquid
const TAU    = ELECTRODE ? 0.2   : 0.5    # transfer relaxation time
const LAM    = 1.0              # theta15, molar volume times csat
const GAM    = 1.0              # sharpening parameter (u = 0, or well above max|u|)
const EPS    = DX               # interface thickness
# Nucleation threshold: J = 0 wherever alpha <= AMIN, so trace gas has no source
# (Example 31). CDI_AMIN=0 is the note's model (the electrode case's default).
const AMIN   = parse(Float64, get(ENV, "CDI_AMIN", ELECTRODE ? "0.0" : "0.01"))

# Seeds. Bubbles: centre x, centre y, radius (A, B, E isolated; C and D the
# close pair, gap 0.06 between the rims). Electrode: site x on the wall, radius.
const BUBBLES = [(0.30, 0.30, 0.08),    # A
                 (0.30, 0.72, 0.06),    # B
                 (0.61, 0.68, 0.06),    # C
                 (0.79, 0.68, 0.06),    # D
                 (0.72, 0.28, 0.04)]    # E
const NAMES = ["A", "B", "C", "D", "E"]
const SITES = [(0.12, 0.05),
               (0.37, 0.07),
               (0.55, 0.04),
               (0.82, 0.06)]
minimum(ELECTRODE ? (s[2] for s in SITES) : (b[3] for b in BUBBLES)) >= 4 * EPS ||
    error("radius < 4 eps: refine the mesh")

const T_END  = parse(Float64, get(ENV, "CDI_T", ELECTRODE ? "2.0" : "1.5"))
# Explicit limits of the note in 2D (dx^2/(4 gamma eps), dx^2/(4 D), tau), with
# a factor for the consistent mass matrix
const NSTEPS = ceil(Int, T_END / (0.2 * min(DX^2 / (4 * GAM * EPS), DX^2 / (4 * D), TAU)) - 1e-9)
const DT     = T_END / NSTEPS   # the step Finch takes (Nsteps is what it is given)

# Every solve is with a mass matrix (or the Laplacian of psi): CG.
const NAME = ELECTRODE ? "cdielectrode" : "cdibubbles"
ctx = init_finch(NAME * "2d"; verbose = false, finch_solver_kwargs(:cg)...,
                 target = finch_target(NAME, PROJ_DIR))

domain(2)
if ELECTRODE
    mesh(MESH_QUAD; elsperdim = NEL, interval = (0.0, L), periodic = :x)
else
    mesh(MESH_QUAD; elsperdim = NEL, interval = (0.0, L))
end
function_space(order = 1)

alpha = variable("alpha")
c     = variable("c")
g     = variable("g"; kind = VAR_VECTOR, components = 2)   # projected gradient of alpha (nodal normal)
J     = variable("J")       # transfer rate (auxiliary field, set by update)
psi   = ELECTRODE ? variable("psi") : nothing   # velocity potential, u = grad(psi)
test_symbol("v")

if ELECTRODE
    boundary(alpha, :ymax, BC_DIRICHLET, 0.0)
    boundary(c,     :ymax, BC_DIRICHLET, CSAT)
    boundary(psi,   :ymax, BC_DIRICHLET, 0.0)
else
    for b in (:xmin, :xmax, :ymin, :ymax)
        boundary(alpha, b, BC_DIRICHLET, 0.0)
        boundary(c,     b, BC_DIRICHLET, CSUPER)
    end
end

# Equilibrium tanh profile around each seed, combined by the maximum
wrap(d) = ifelse(d > L/2, d - L, ifelse(d < -L/2, d + L, d))   # traced: no ?:
seed_bubbles(x, y) = maximum(0.5 * (1 - tanh((sqrt((x - bx)^2 + (y - by)^2) - r) / (2 * EPS)))
                             for (bx, by, r) in BUBBLES)
seed_sites(x, y)   = maximum(0.5 * (1 - tanh((sqrt(wrap(x - sx)^2 + y^2) - r) / (2 * EPS)))
                             for (sx, r) in SITES)
seed(x, y) = ELECTRODE ? seed_sites(x, y) : seed_bubbles(x, y)
const C0 = ELECTRODE ? CSAT : CSUPER      # the liquid: saturated / supersaturated
initial(alpha, (x, y, z, t) -> seed(x, y))
initial(c,     (x, y, z, t) -> C0 * (1 - seed(x, y)))

coefficient("gam",   GAM)
coefficient("eps",   EPS)
coefficient("D",     D)
coefficient("Deps",  D / EPS)   # the form language has no division
coefficient("Lam",   LAM)
coefficient("q",     Q)
coefficient("delta", 1e-12)

time_interval(T_END)
# bubbles: ~50 frames per unit time; electrode: ~100 frames
const SNAP_EVERY = ELECTRODE ? max(1, NSTEPS ÷ 100) : max(1, round(Int, NSTEPS / (50 * T_END)))
time_stepper(STEP_EULER_EXPLICIT; Nsteps = NSTEPS, store_trajectory = true, snapshot_every = SNAP_EVERY)

# The expansion flow carries both fields (conservative form); the electrode
# flux q enters the c equation through the bottom wall
const ADV_C     = ELECTRODE ? " - c*dot(grad(psi), grad(v))" : ""
const ADV_ALPHA = ELECTRODE ? " - alpha*dot(grad(psi), grad(v))" : ""
const FLUX_C    = ELECTRODE ? " - q*v*ds(:ymin)" : ""

ELECTRODE && weak_form(psi, "dot(grad(psi), grad(v)) + Lam*J*v")
weak_form(g,   "dot(g, v) - dot(grad(alpha), v)")
weak_form(c,   "Dt(c)*v" * ADV_C * " + D*dot(grad(c), grad(v))" *
               " + Deps*alpha*c*" * NDOTGV_NODAL_2D * " + J*v" * FLUX_C)
weak_form(alpha, "Dt(alpha)*v" * ADV_ALPHA *
                 " + gam*eps*dot(grad(alpha), grad(v))" *
                 " - gam*(alpha - alpha*alpha)*" * NDOTGV_NODAL_2D * " - Lam*J*v")

# Gas volume of each bubble: alpha integrated over its own region of the domain
region(name) = Dict("A" => "ifelse(x < 0.5, ifelse(y < 0.5, alpha, 0), 0)",
                    "B" => "ifelse(x < 0.5, ifelse(y > 0.5, alpha, 0), 0)",
                    "C" => "ifelse(x > 0.5, ifelse(x < 0.7, ifelse(y > 0.5, alpha, 0), 0), 0)",
                    "D" => "ifelse(x > 0.7, ifelse(y > 0.5, alpha, 0), 0)",
                    "E" => "ifelse(x > 0.5, ifelse(y < 0.5, alpha, 0), 0)")[name]

entries = Any[update(J => "ifelse(alpha > $AMIN, ifelse(c - $CSAT*(1 - alpha) > 0, alpha*$(1/TAU)*(c - $CSAT*(1 - alpha)), 0), 0)")]
ELECTRODE && push!(entries, psi)
append!(entries, [g, c, alpha,
                  monitor("gas_volume", "alpha"),
                  monitor("created_rate", "Lam*J"),
                  monitor("gas_moles", "c + alpha*$(1/LAM)")])
ELECTRODE || append!(entries, [monitor("vol_" * n, region(n)) for n in NAMES])
step_schedule(entries...)

println("NEL = $NEL, dt = $DT, steps = $NSTEPS, AMIN = $AMIN")
t0 = time()
solve()
ELECTRODE ? finch_run!(ctx, alpha, c, psi) : finch_run!(ctx, alpha, c)
println("run time = ", round(time() - t0, digits = 1), " s")

# Snapshots from Finch's results: columns 1-2 are x, y, then one per snapshot
x = ctx.grid.nodes[1, :]; y = ctx.grid.nodes[2, :]
results = joinpath(PROJ_DIR, "results")
times, fields = read_snapshots(PROJ_DIR)
check_final_time(times, T_END)
A = fields["alpha"]
C = fields["c"]
write_timeseries(OUT_DIR, "alpha", [x y], A)
write_timeseries(OUT_DIR, "c", [x y], C)
ELECTRODE && write_timeseries(OUT_DIR, "psi", [x y], fields["psi"])
write_snapshot_times(OUT_DIR, times)

# history.csv, every step: t, int alpha, volume created so far, int c + int
# alpha/theta15, then the gas volume of bubbles A..E (bubbles) or the gas the
# electrode has put in, q*L*t (electrode). Row k of created_rate is the J step
# k used (row 0, before the first update, is zero).
mon = readdlm(joinpath(results, "monitors.csv"), ',', skipstart = 1)
t_m, vol, rate, moles = mon[:, 2], mon[:, 3], mon[:, 4], mon[:, 5]
created = DT .* cumsum([0.0; rate[2:end]])
writedlm(joinpath(OUT_DIR, "history.csv"),
         ELECTRODE ? [t_m vol created moles Q*L .* t_m] : [t_m vol created moles mon[:, 6:end]], ',')

a, cc = A[:, end], C[:, end]
println("gas volume: ", vol[1], " -> ", vol[end], ", created ", created[end],
        ", balance error ", (vol[end] - vol[1]) - created[end])
if ELECTRODE
    far = maximum(a[k] for k in eachindex(x) if y[k] > 0.6)
    println("gas moles (dissolved + gas): ", moles[1], " -> ", moles[end],
            ", electrode put in ", Q*L*T_END, " (the rest left through the top)")
else
    dist(k) = minimum(sqrt((x[k] - bx)^2 + (y[k] - by)^2) - r for (bx, by, r) in BUBBLES)
    far = maximum(a[k] for k in eachindex(x) if dist(k) > 0.15)
    println("gas moles (dissolved + gas): ", moles[1], " -> ", moles[end],
            ", change ", moles[end] - moles[1], " (diffusion through the boundary)")
    for (i, n) in enumerate(NAMES)
        v0, v1 = mon[1, 5 + i], mon[end, 5 + i]
        println("bubble $n: equivalent radius ", round(sqrt(v0 / pi), digits = 4), " -> ",
                round(sqrt(v1 / pi), digits = 4))
    end
end
println("min alpha = ", minimum(a), ", max alpha = ", maximum(a))
println("min c = ", minimum(cc), ", max c = ", maximum(cc))
println(ELECTRODE ? "max alpha in the upper liquid (y > 0.6) = " : "max alpha away from the bubbles = ", far)
println(length(times), " snapshots written")
