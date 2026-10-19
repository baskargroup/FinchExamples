#==============================================================================#
# Sharpening-term study for the semi-implicit 1D CDI scheme of main.jl.
#
# The sharpening flux gamma*(phi - phi^2)*nhat is split between the unknown phi
# and a lagged field phis in several ways. All of them fit one weak form,
#     gamma*[(c1 + c2*phis)*phi + (c3 + c4*phis)*phis]*nhat(phis),
# whose four numbers are Julia globals read at assembly (like the BDF
# coefficients), so one Finch setup runs every form:
#   lagged    phis*(1 - phis)              c = (0, 0, 1, -1)  right-hand side only
#   derived   phi*(1 - phis)               c = (1, -1, 0, 0)  as in main.jl
#   reversed  phis*(1 - phi)               c = (0, -1, 1, 0)
#   newton    (1 - 2*phis)*phi + phis^2    c = (1, -2, 0, 1)  Newton about phis
#   extrap    pe*(1 - pe), pe = 2*phio - phioo, one solve per step
#             (BDF2 with an extrapolated explicit term, the usual IMEX choice)
# The first four run with 1, 2 or 3 Picard iterations: iteration 1 lags with
# phio, each later one with the previous iterate. 'lagged' with one iteration
# is the explicit treatment phio*(1 - phio). The normal is always recovered
# from phis (see main.jl).
#
# Each form runs at six CFL numbers on the mesh given on the command line, with
# eps = dx. dt is rounded down so that T is reached exactly. Measure: the max
# nodal error against the exact profile at T. A run that produces NaN or
# leaves [-1, 2] is stopped and marked diverged.
#
# Writes into output/sharpening_sweep/ (runs already in results.csv are skipped):
#   results.csv          run, nel, form, iterations, CFL, dt, steps, status,
#                        t_end, max error, min phi, max phi, seconds
#   hist/<run>.csv       step, t, max error, min phi, max phi
#   frames/<run>.bin     phi at 41 instants for the movies: Int64 nnode, Int64 nf,
#                        Float64 times[nf], Float32 phi[nnode, nf]
#   nodes_n<nel>.csv     node coordinates
#
# Usage:  julia sharpening_sweep.jl nel [forms] [iterations] [CFLs]
#   e.g.  julia sharpening_sweep.jl 200 derived 2 0.25
#==============================================================================#

using Finch, DelimitedFiles, Printf

nel = parse(Int, ARGS[1])
listarg(i, f, default) = length(ARGS) >= i ? f.(split(ARGS[i], ",")) : default
FORMS = listarg(2, String, ["extrap", "lagged", "derived", "reversed", "newton"])
ITERS = listarg(3, s -> parse(Int, s), [1, 2, 3])
CFLS  = listarg(4, s -> parse(Float64, s), [2.0, 1.0, 0.5, 0.25, 0.125, 0.0625])

outdir = "output/sharpening_sweep"
mkpath(joinpath(outdir, "hist")); mkpath(joinpath(outdir, "frames"))

initFinch("cdi1dsweep");
useLog("sweeplog_n$nel", dir=outdir, level=1)

domain(1)
functionSpace(order=1)
mesh(LINEMESH, elsperdim=nel, bids=2)

phi   = variable("phi")
phio  = variable("phio")
phioo = variable("phioo")
phis  = variable("phis")
nhat  = variable("nhat")
testSymbol("v")

# Same case as main.jl
dx   = 1.0/nel
uvel = 1.0
gam  = 1.0
eps  = dx
x0   = 0.25
T    = 0.25

boundary(phi, 1, DIRICHLET, 0)
boundary(phi, 2, DIRICHLET, 1)
initial(phi, "0.5*(1 + tanh((x - $x0)/(2*$eps)))")

# BDF and sharpening-form coefficients, read by Finch at every assembly
a0_dt = 0.0;  a1_dt = 0.0;  a2_dt = 0.0
sh_c1 = 0.0;  sh_c2 = 0.0;  sh_c3 = 0.0;  sh_c4 = 0.0
const FORM_C = Dict("lagged"   => (0.0,  0.0, 1.0, -1.0),
                    "extrap"   => (0.0,  0.0, 1.0, -1.0),
                    "derived"  => (1.0, -1.0, 0.0,  0.0),
                    "reversed" => (0.0, -1.0, 1.0,  0.0),
                    "newton"   => (1.0, -2.0, 0.0,  1.0))

coefficient("theta1", 1.0)
coefficient("theta2", 1.0)
coefficient("theta3", 1.0)
coefficient("theta4", 1.0)
coefficient("u", [string(uvel)], type=VECTOR)
coefficient("gamma", gam)
coefficient("epsilon", eps)
coefficient("a0dt", "Main.a0_dt")
coefficient("a1dt", "Main.a1_dt")
coefficient("a2dt", "Main.a2_dt")
coefficient("c1", "Main.sh_c1")
coefficient("c2", "Main.sh_c2")
coefficient("c3", "Main.sh_c3")
coefficient("c4", "Main.sh_c4")

# Nodal normal by gradient recovery, identical to main.jl
const g_grid  = Finch.finch_state.grid_data
const g_refel = Finch.finch_state.refel
const g_geo   = Finch.finch_state.geo_factors
const g_Np    = g_refel.Np
const g_Nqp   = g_refel.Nqp
const g_nnode = size(g_grid.allnodes, 2)
const g_nel   = size(g_grid.loc2glb, 2)
const g_lump  = zeros(g_nnode)
for e = 1:g_nel
    detj = g_geo.detJ[1, e]
    for i = 1:g_Np
        acc = 0.0
        for q = 1:g_Nqp
            acc += g_refel.wg[q] * g_refel.Q[q, i]
        end
        g_lump[g_grid.loc2glb[i, e]] += acc * detj
    end
end
const g_gx = zeros(g_nnode)
const g_pe = zeros(g_Np)
const g_dr = zeros(g_Nqp)

function recover_nhat!(nhatvals, p, delta)
    fill!(g_gx, 0.0)
    @inbounds for e = 1:g_nel
        rx = g_geo.J[e].rx[1]
        detj = g_geo.detJ[1, e]
        for i = 1:g_Np
            g_pe[i] = p[g_grid.loc2glb[i, e]]
        end
        for q = 1:g_Nqp
            dr = 0.0
            for i = 1:g_Np
                dr += g_refel.Qr[q, i] * g_pe[i]
            end
            g_dr[q] = dr
        end
        for i = 1:g_Np
            ax = 0.0
            for q = 1:g_Nqp
                ax += g_refel.wg[q] * g_refel.Q[q, i] * g_dr[q] * rx
            end
            g_gx[g_grid.loc2glb[i, e]] += ax * detj
        end
    end
    @inbounds for k = 1:g_nnode
        gx = g_gx[k] / g_lump[k]
        nhatvals[k] = gx / sqrt(gx*gx + delta)
    end
    return nothing
end

# main.jl's weak form with the sharpening term in the general form above
weakForm(phi, "theta1*a0dt*phi*v" *
              " + theta1*(a1dt*phio + a2dt*phioo)*v" *
              " - 0.5*theta2*phi*dot(u, grad(v))" *
              " - 0.5*theta2*phi*(grad(u)*v + dot(u, grad(v)))" *
              " + theta3*gamma*epsilon*dot(grad(phi), grad(v))" *
              " - theta4*gamma*(c1 + c2*phis)*phi .* nhat*grad(v)" *
              " - theta4*gamma*(c3 + c4*phis)*phis .* nhat*grad(v)")

evalInitialConditions()
phi0 = copy(phi.values)
x = g_grid.allnodes[1, :]
exact(t) = 0.5 .* (1 .+ tanh.((x .- x0 .- uvel*t) ./ (2*eps)))
writedlm(joinpath(outdir, "nodes_n$nel.csv"), x, ',')

function set_bdf!(dt, order)
    if order == 1
        global a0_dt = 1.0/dt; global a1_dt = -1.0/dt; global a2_dt = 0.0
    else
        global a0_dt = 1.5/dt; global a1_dt = -2.0/dt; global a2_dt = 0.5/dt
    end
end
function set_form!(form)
    global sh_c1, sh_c2, sh_c3, sh_c4 = FORM_C[form]
end

function solve_with_lag!()
    recover_nhat!(nhat.values, phis.values, 1e-12)
    solve(phi)
end

const NF = 40      # movie frames after t = 0

# One run to t = T. Returns the results row; writes the history and frames.
function run_case(form, iters, CFL, name)
    nsteps = ceil(Int, T/(CFL*dx/uvel) - 1e-9)
    dt = T/nsteps
    set_form!(form)
    set_bdf!(dt, 1)
    phi.values .= phi0; phio.values .= phi0; phioo.values .= phi0; phis.values .= phi0
    fsteps = [round(Int, k*nsteps/NF) for k = 0:NF]
    frames = Matrix{Float32}(undef, g_nnode, NF+1); ftimes = Float64[]
    function store_frames!(step)
        for k in findall(==(step), fsteps)
            push!(ftimes, step*dt); frames[:, length(ftimes)] .= phi.values[:]
        end
    end
    store_frames!(0)
    hist = zeros(0, 5)
    status = "ok"; t_end = 0.0; err = NaN
    t0 = time()
    for step = 1:nsteps
        step == 2 && set_bdf!(dt, 2)
        if form == "extrap"
            if step == 1
                phis.values .= phio.values
            else
                phis.values .= 2 .* phio.values .- phioo.values
            end
            solve_with_lag!()
        else
            phis.values .= phio.values
            solve_with_lag!()
            for it = 2:iters
                phis.values .= phi.values
                solve_with_lag!()
            end
        end
        phioo.values .= phio.values
        phio.values  .= phi.values
        t_end = step*dt
        p = phi.values[:]
        if !all(isfinite, p) || maximum(p) > 2 || minimum(p) < -1
            status = "diverged"
            break
        end
        err = maximum(abs.(p .- exact(t_end)))
        hist = vcat(hist, [step t_end err minimum(p) maximum(p)])
        store_frames!(step)
    end
    secs = time() - t0
    p = phi.values[:]
    if status == "ok" && (minimum(p) < -0.01 || maximum(p) > 1.01)
        status = "overshoot"            # finished, but phi left [-0.01, 1.01]
    end
    status == "diverged" && (err = NaN)
    writedlm(joinpath(outdir, "hist", name * ".csv"), hist, ',')
    open(joinpath(outdir, "frames", name * ".bin"), "w") do io
        nf = length(ftimes)
        write(io, Int64(g_nnode), Int64(nf), ftimes, frames[:, 1:nf])
    end
    return [name, nel, form, iters, CFL, dt, nsteps, status, t_end, err,
            minimum(p), maximum(p), round(secs, digits=2)]
end

resfile = joinpath(outdir, "results.csv")
if !isfile(resfile)
    open(resfile, "w") do io
        println(io, "run,nel,form,iterations,CFL,dt,steps,status,t_end,maxerr,phimin,phimax,seconds")
    end
end
done = Set(String(split(l, ",")[1]) for l in readlines(resfile)[2:end])

for CFL in CFLS, form in FORMS, iters in ITERS
    form == "extrap" && iters != 1 && continue
    name = @sprintf("n%d_%s_it%d_cfl%g", nel, form, iters, CFL)
    name in done && continue
    row = run_case(form, iters, CFL, name)
    open(resfile, "a") do io
        println(io, join(row, ","))
    end
    @printf("%-28s %-9s t_end %.4f  max error %.4e  phi [%.4f, %.4f]  %.1f s\n",
            name, row[8], row[9], row[10], row[11], row[12], row[13])
    flush(stdout)
end

finalizeFinch()
