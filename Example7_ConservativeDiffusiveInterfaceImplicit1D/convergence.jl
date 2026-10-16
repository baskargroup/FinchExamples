#==============================================================================#
# Temporal convergence study for the semi-implicit 1D CDI scheme of main.jl.
#
# The error against the exact tanh contains a spatial part (~5e-3 on this mesh)
# that does not shrink with dt, so it cannot show the order in time. Here each
# run is instead compared with a reference solution on the SAME mesh at a very
# small dt: the spatial error is identical in both and cancels, and what is
# left is the time-stepping error alone.
#
# Two variants are measured: the scheme as derived (two Picard iterations per
# step, expected second order) and a single-iteration variant (lagged field =
# previous level, expected first order).
#
# Same case and discretization as main.jl. Writes output/convergence.csv with
# columns  CFL, dt, error(2 iterations), error(1 iteration).
#==============================================================================#

using Finch

initFinch("cdi1dconv");
useLog("cdi1dconvlog", dir="output", level=1)

domain(1)
functionSpace(order=1)

nel = 200
mesh(LINEMESH, elsperdim=nel, bids=2)

phi   = variable("phi")
phio  = variable("phio")
phioo = variable("phioo")
phis  = variable("phis")
testSymbol("v")

dx   = 1.0/nel
uvel = 1.0
gam  = 1.0
eps  = dx
x0   = 0.25
T    = 0.25

boundary(phi, 1, DIRICHLET, 0)
boundary(phi, 2, DIRICHLET, 1)
initial(phi, "0.5*(1 + tanh((x - $x0)/(2*$eps)))")

# BDF coefficients over dt, read by Finch at every assembly (see main.jl)
a0_dt = 0.0;  a1_dt = 0.0;  a2_dt = 0.0

coefficient("theta1", 1.0)
coefficient("theta2", 1.0)
coefficient("theta3", 1.0)
coefficient("theta4", 1.0)
coefficient("u", [string(uvel)], type=VECTOR)
coefficient("gamma", gam)
coefficient("epsilon", eps)
coefficient("delta", 1e-12)
coefficient("a0dt", "Main.a0_dt")
coefficient("a1dt", "Main.a1_dt")
coefficient("a2dt", "Main.a2_dt")

parameter("nhat", "grad(phis) ./ sqrt.(dot(grad(phis), grad(phis)) .+ delta)")

# Identical to main.jl
weakForm(phi, "theta1*a0dt*phi*v" *
              " + theta1*(a1dt*phio + a2dt*phioo)*v" *
              " - 0.5*theta2*phi*dot(u, grad(v))" *
              " - 0.5*theta2*phi*(grad(u)*v + dot(u, grad(v)))" *
              " + theta3*gamma*epsilon*dot(grad(phi), grad(v))" *
              " - theta4*gamma*phi*(1 - phis) .* dot(nhat, grad(v))")

evalInitialConditions()
phi0 = copy(phi.values)      # the initial field, restored before every run

# Run the scheme to t = T with the given CFL and number of Picard iterations
# (2 = as derived, 1 = lagged field is the previous level). Returns phi(T).
function run_case(CFL, iters)
    dt = CFL*dx/uvel
    nsteps = Int(round(T/dt))
    @assert abs(nsteps*dt - T) < 1e-12 "T/dt must be an integer"
    phi.values .= phi0; phio.values .= phi0; phioo.values .= phi0; phis.values .= phi0
    global a0_dt = 1.0/dt; global a1_dt = -1.0/dt; global a2_dt = 0.0     # BDF1
    for step = 1:nsteps
        if step == 2
            global a0_dt = 1.5/dt; global a1_dt = -2.0/dt; global a2_dt = 0.5/dt   # BDF2
        end
        phis.values .= phio.values
        solve(phi)
        if iters == 2
            phis.values .= phi.values
            solve(phi)
        end
        phioo.values .= phio.values
        phio.values  .= phi.values
    end
    return copy(phi.values[:])
end

# Reference: the derived scheme at a dt far below the sweep
CFL_ref = 1/64
println("reference run at CFL $CFL_ref ...")
phi_ref = run_case(CFL_ref, 2)

CFLs = [1.0, 0.5, 0.25, 0.125, 0.0625]
err2 = Float64[]; err1 = Float64[]
for c in CFLs
    push!(err2, maximum(abs.(run_case(c, 2) .- phi_ref)))
    push!(err1, maximum(abs.(run_case(c, 1) .- phi_ref)))
    println("CFL $(c):  dt = $(c*dx)   2 iterations: $(err2[end])   1 iteration: $(err1[end])")
end

finalizeFinch()

# Observed order from a least-squares fit of log(error) against log(dt)
slope(e) = (lx = log.(CFLs .* dx); ly = log.(e); n = length(lx);
            (n*sum(lx .* ly) - sum(lx)*sum(ly)) / (n*sum(lx .^ 2) - sum(lx)^2))
println("fitted order, 2 iterations: ", round(slope(err2), digits=2))
println("fitted order, 1 iteration:  ", round(slope(err1), digits=2))

using DelimitedFiles
writedlm("output/convergence.csv", [CFLs CFLs .* dx err2 err1], ',')
println("output/convergence.csv written (CFL, dt, err 2 iter, err 1 iter)")
