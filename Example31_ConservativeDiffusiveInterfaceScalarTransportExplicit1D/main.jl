#==============================================================================#
# 1D CDI phase field with dissolved-gas transport: growth of a gas column,
# solved explicitly with Finch V2 (budvinchathura/Finch-V2, branch
# mlir-batched-linalg), Julia target. The formulation is the P30 note "CDI with
# Scalar Transport, Explicit Galerkin FE Formulation for Diffusive Bubble
# Growth", after Qiu and Mani (2026):
#
#   alpha_t + (u*alpha)_x - (gamma*eps*alpha_x)_x
#     + (gamma*(alpha - alpha^2)*nhat)_x - theta15*J = 0
#
#   c_t + (u*c)_x - (D*c_x)_x - (D*alpha*c/eps*nhat)_x + J = 0
#
#   J = alpha/tau * max(0, c - csat*(1 - alpha)),     u_x = theta15*J
#
# alpha is the gas volume fraction (1 in the gas), c the dissolved gas per unit
# total volume, nhat = alpha_x/|alpha_x| the interface normal. Gas crossing the
# interface at the rate J expands the column; theta15 = R*T/P is the molar
# volume. The fourth term of the c equation keeps the dissolved gas out of the gas.
#
# One forward Euler step is a step_schedule, every entry reading the previous level:
#   1. update J    node by node
#   2. solve  psi  : psi_xx = theta15*J, psi = 0 at both ends      quasi-static
#      The velocity is u = psi_x (no momentum equation). J is symmetric about
#      the centre and psi = 0 at both ends, so u(L/2) = 0: the column stays in
#      place and grows symmetrically.
#   3. solve  gx   : (gx, v) = (alpha_x, v)                         quasi-static
#      the nodal gradient by L2 projection (consistent mass),
#      nhat = gx/sqrt(gx^2 + delta)
#   4. step   c    (reads alpha, psi, gx and J of the previous level)
#   5. step   alpha
# c goes before alpha so that it reads the previous alpha, not the new one.
# Advection is in conservative form: u is not divergence free here, and the
# conservative form keeps d/dt int(alpha) = int(theta15*J) exact.
#
# Test problem (Qiu and Mani's 1D demonstration): a column of width 1/4 at the
# centre of [0, 1], liquid supersaturated at c = 3*csat, alpha = 0 and
# c = csuper held at both ends. The slides give only csuper and the width; D,
# tau, gamma, eps and the molar volume below are this example's choice.
# Non-dimensional: L = 1, c scaled by csat, theta15 = Lambda = R*T*csat/P.
#
# Trace gas in the supersaturated liquid grows like exp(2*theta15*t/tau), since
# J is proportional to alpha; that includes the tanh tails, which are never
# exactly zero. With theta15/tau = 10 and T = 1 small gas pockets appear ahead
# of the interfaces (alpha ~ 3e-3); at theta15/tau = 20 gas fills the liquid.
# CDI_AMIN sets a nucleation threshold: J = 0 wherever alpha <= AMIN, so the
# trace gas has no source and stays at its initial level. The default 0.01
# removes the pockets (max alpha ahead of the column 3e-3 -> 1e-9) and changes
# the column edges by 6e-4; CDI_AMIN=0 is the note's model.
#
# Continuous Galerkin, linear elements. Run with ./run.sh (Julia 1.12).
#==============================================================================#

include(joinpath(@__DIR__, "..", "common.jl"))   # Finch, FINCH_TARGET switch, helpers

# CSV results go to output/; Finch's generated project (and its raw results)
# to output/finch/
const OUT_DIR  = joinpath(@__DIR__, "output")
const PROJ_DIR = joinpath(OUT_DIR, "finch")
mkpath(PROJ_DIR)

# Parameters
const NEL    = parse(Int, get(ENV, "CDI_NEL", "200"))
const L      = 1.0
const AMIN   = parse(Float64, get(ENV, "CDI_AMIN", "0.01"))   # nucleation threshold, see above
const DX     = L / NEL
const CSAT   = 1.0
const CSUPER = 3.0 * CSAT
const W      = 0.25 * L
const D      = 0.05
const TAU    = 0.1
const LAM    = 1.0              # theta15, molar volume times csat
const GAM    = 1.0
const EPS    = DX

const T_END  = 1.0
const NSTEPS = round(Int, T_END / (0.2 * min(DX^2 / (2 * GAM * EPS), DX^2 / (2 * D), TAU)))
const DT     = T_END / NSTEPS   # the step Finch takes (Nsteps is what it is given)

# Every solve here is SPD (mass matrices and a Laplacian): CG.
ctx = init_finch("cdigas1d"; verbose = false, finch_solver_kwargs(:cg)...,
                 target = finch_target("cdigas", PROJ_DIR))

domain(1)
mesh(MESH_LINE; elsperdim = NEL, interval = (0.0, L))
function_space(order = 1)

alpha = variable("alpha")
c     = variable("c")
psi   = variable("psi")     # velocity potential, u = grad(psi)
gx    = variable("gx")      # projected gradient of alpha
J     = variable("J")       # transfer rate (auxiliary field, set by update)
test_symbol("v")

for b in (:xmin, :xmax)
    boundary(alpha, b, BC_DIRICHLET, 0.0)
    boundary(c,     b, BC_DIRICHLET, CSUPER)
    boundary(psi,   b, BC_DIRICHLET, 0.0)
end

column(x) = 0.5 * (tanh((x - (L - W) / 2) / (2 * EPS)) - tanh((x - (L + W) / 2) / (2 * EPS)))
initial(alpha, (x, y, z, t) -> column(x))
initial(c,     (x, y, z, t) -> CSUPER * (1 - column(x)))

coefficient("gam",   GAM)
coefficient("eps",   EPS)
coefficient("D",     D)
coefficient("Deps",  D / EPS)   # the form language has no division
coefficient("Lam",   LAM)
coefficient("delta", 1e-12)

time_interval(T_END)
time_stepper(STEP_EULER_EXPLICIT; Nsteps = NSTEPS, store_trajectory = true,
             snapshot_every = NSTEPS ÷ 20)

weak_form(psi, "dot(grad(psi), grad(v)) + Lam*J*v")
weak_form(gx,  "gx*v - grad(alpha)*v")
weak_form(c,   "Dt(c)*v - c*dot(grad(psi), grad(v)) + D*dot(grad(c), grad(v))" *
               " + Deps*alpha*c*" * NDOTGV_NODAL_1D * "*grad(v) + J*v")
weak_form(alpha, "Dt(alpha)*v - alpha*dot(grad(psi), grad(v))" *
                 " + gam*eps*dot(grad(alpha), grad(v))" *
                 " - gam*(alpha - alpha*alpha)*" * NDOTGV_NODAL_1D * "*grad(v)" *
                 " - Lam*J*v")

step_schedule(update(J => "ifelse(alpha > $AMIN, ifelse(c - $CSAT*(1 - alpha) > 0, alpha*$(1/TAU)*(c - $CSAT*(1 - alpha)), 0), 0)"),
              psi, gx, c, alpha,
              monitor("gas_volume", "alpha"),
              monitor("created_rate", "Lam*J"),
              monitor("gas_moles", "c + alpha*$(1/LAM)"))

println("dt = $DT, steps = $NSTEPS, AMIN = $AMIN")
t0 = time()
solve()
finch_run!(ctx, alpha, c, psi)
println("run time = ", round(time() - t0, digits = 1), " s")

# Nodes in order of x
x  = ctx.grid.nodes[1, :]
o  = sortperm(x)
xs = x[o]

# u = psi_x: constant on each element; at a node, the average of its two elements
function velocity(ps)
    slope = diff(ps) ./ diff(xs)
    return [slope[1]; 0.5 .* (slope[1:end-1] .+ slope[2:end]); slope[end]]
end

# alpha = 0.5 crossings, left and right edges of the column
function edges(a)
    i = findfirst(>=(0.5), a); j = findlast(>=(0.5), a)
    xl = xs[i-1] + (0.5 - a[i-1]) / (a[i] - a[i-1]) * (xs[i] - xs[i-1])
    xr = xs[j] + (0.5 - a[j]) / (a[j+1] - a[j]) * (xs[j+1] - xs[j])
    return xl, xr
end

# Snapshots from Finch's results: one column per snapshot, x first
results = joinpath(PROJ_DIR, "results")
times, fields = read_snapshots(PROJ_DIR)
check_final_time(times, T_END)
A = fields["alpha"][o, :]
C = fields["c"][o, :]
U = reduce(hcat, [velocity(fields["psi"][o, k]) for k in eachindex(times)])
write_timeseries(OUT_DIR, "alpha", reshape(xs, :, 1), A)
write_timeseries(OUT_DIR, "c", reshape(xs, :, 1), C)
write_timeseries(OUT_DIR, "u", reshape(xs, :, 1), U)
write_snapshot_times(OUT_DIR, times)
E = [edges(A[:, k]) for k in axes(A, 2)]
writedlm(joinpath(OUT_DIR, "edges.csv"), [times first.(E) last.(E)], ',')

a, cc, u = A[:, end], C[:, end], U[:, end]
writedlm(joinpath(OUT_DIR, "solution.csv"), [xs a cc u], ',')

# Gas volume balance from the monitors: int alpha(t) - int alpha(0) against
# the sum of dt * int theta15*J. Row k of created_rate is the J step k used
# (row 0, before the first update, is zero).
# history.csv: t, int alpha, volume created so far, int c + int alpha/theta15
mon = readdlm(joinpath(results, "monitors.csv"), ',', skipstart = 1)
t_m, vol, rate, moles = mon[:, 2], mon[:, 3], mon[:, 4], mon[:, 5]
created = DT .* cumsum([0.0; rate[2:end]])
writedlm(joinpath(OUT_DIR, "history.csv"), [t_m vol created moles], ',')
far = maximum(A[k, j] for k in eachindex(xs), j in axes(A, 2) if abs(xs[k] - L/2) > 0.4)

xl, xr = E[end]
println("column: [", round(E[1][1], digits = 4), ", ", round(E[1][2], digits = 4), "] -> [",
        round(xl, digits = 4), ", ", round(xr, digits = 4), "]")
println("gas volume: ", vol[1], " -> ", vol[end], ", created ", created[end],
        ", balance error ", (vol[end] - vol[1]) - created[end])
println("max |u| = ", maximum(abs.(u)), " (gamma = $GAM)")
println("min alpha = ", minimum(a), ", max alpha = ", maximum(a))
println("min c = ", minimum(cc), ", max c = ", maximum(cc))
println("max alpha away from the column = ", far)
println(length(times), " snapshots written")
