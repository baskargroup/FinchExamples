#==============================================================================#
# Sharpening-term study for the semi-implicit 2D CDI scheme of main.jl.
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
# from phis, and the outflow flux uses phis, as in main.jl.
#
# Each form runs at six CFL numbers on the mesh given on the command line, with
# eps = dx. dt is rounded down so that T is reached exactly. Measure: the max
# nodal error against the translated tanh circle at T (not an exact solution
# for a curved interface, so the error has a floor set by the mesh). A run that
# produces NaN or leaves [-1, 2] is stopped and marked diverged.
#
# Writes into output/sharpening_sweep/ (runs already in results.csv are skipped):
#   results.csv          run, nel, form, iterations, CFL, dt, steps, status,
#                        t_end, max error, min phi, max phi, relative mass
#                        change, seconds
#   hist/<run>.csv       step, t, max error, min phi, max phi, mass change
#   frames/<run>.bin     phi at 41 instants for the movies: Int64 nnode, Int64 nf,
#                        Float64 times[nf], Float32 phi[nnode, nf]
#   nodes_n<nel>.csv     node coordinates x, y
#
# Usage:  julia sharpening_sweep.jl nel [forms] [iterations] [CFLs]
#   e.g.  julia sharpening_sweep.jl 100 derived 2 0.25
#==============================================================================#

using Finch, DelimitedFiles, Printf

nel = parse(Int, ARGS[1])
listarg(i, f, default) = length(ARGS) >= i ? f.(split(ARGS[i], ",")) : default
FORMS = listarg(2, String, ["extrap", "lagged", "derived", "reversed", "newton"])
ITERS = listarg(3, s -> parse(Int, s), [1, 2, 3])
CFLS  = listarg(4, s -> parse(Float64, s), [2.0, 1.0, 0.5, 0.25, 0.125, 0.0625])

outdir = "output/sharpening_sweep"
mkpath(joinpath(outdir, "hist")); mkpath(joinpath(outdir, "frames"))

initFinch("cdi2dsweep");
useLog("sweeplog_n$nel", dir=outdir, level=1)

domain(2)
functionSpace(order=1)
linAlgOptions(iterative=true, method="GMRES", pc="NONE", reltol=1e-12)

mesh(QUADMESH, elsperdim=nel)
addBoundaryID(1, "x==0")           # left,   inflow
addBoundaryID(2, "x==1")           # right,  outflow
addBoundaryID(3, "y==0")
addBoundaryID(4, "y==1")

phi   = variable("phi")
phio  = variable("phio")
phioo = variable("phioo")
phis  = variable("phis")
nhat  = variable("nhat", type=VECTOR)
testSymbol("v")

# Same case as main.jl
dx   = 1.0/nel
ux   = 1.0
uy   = 0.0
umax = sqrt(ux^2 + uy^2)
gam  = umax
eps  = dx
R    = 0.15
xc0  = 0.3
yc0  = 0.5
T    = 0.4

boundary(phi, 1, DIRICHLET, 0)
initial(phi, "0.5*(1 - tanh((sqrt((x-$xc0)^2 + (y-$yc0)^2) - $R)/(2*$eps)))")

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
coefficient("u", [string(ux), string(uy)], type=VECTOR)
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
const g_lump = zeros(g_nnode)
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
const g_gy = zeros(g_nnode)
const g_pe = zeros(g_Np)
const g_dr = zeros(g_Nqp)
const g_ds = zeros(g_Nqp)

function recover_nhat!(nhatvals, p, delta)
    fill!(g_gx, 0.0); fill!(g_gy, 0.0)
    @inbounds for e = 1:g_nel
        J = g_geo.J[e]
        rx = J.rx[1]; ry = J.ry[1]; sx = J.sx[1]; sy = J.sy[1]
        detj = g_geo.detJ[1, e]
        for i = 1:g_Np
            g_pe[i] = p[g_grid.loc2glb[i, e]]
        end
        for q = 1:g_Nqp
            dr = 0.0; ds = 0.0
            for i = 1:g_Np
                dr += g_refel.Qr[q, i] * g_pe[i]
                ds += g_refel.Qs[q, i] * g_pe[i]
            end
            g_dr[q] = dr; g_ds[q] = ds
        end
        for i = 1:g_Np
            ax = 0.0; ay = 0.0
            for q = 1:g_Nqp
                w  = g_refel.wg[q] * g_refel.Q[q, i]
                ax += w * (g_dr[q]*rx + g_ds[q]*sx)
                ay += w * (g_dr[q]*ry + g_ds[q]*sy)
            end
            k = g_grid.loc2glb[i, e]
            g_gx[k] += ax * detj
            g_gy[k] += ay * detj
        end
    end
    @inbounds for k = 1:g_nnode
        gx = g_gx[k] / g_lump[k]
        gy = g_gy[k] / g_lump[k]
        m  = sqrt(gx*gx + gy*gy + delta)
        nhatvals[1, k] = gx/m
        nhatvals[2, k] = gy/m
    end
    return nothing
end

parameter("g", "theta2*phis*dot(u, normal())")

# main.jl's weak form with the sharpening term in the general form above
weakForm(phi, "theta1*a0dt*phi*v" *
              " + theta1*(a1dt*phio + a2dt*phioo)*v" *
              " - 0.5*theta2*phi*dot(u, grad(v))" *
              " - 0.5*theta2*phi*(div(u)*v + dot(u, grad(v)))" *
              " + theta3*gamma*epsilon*dot(grad(phi), grad(v))" *
              " - theta4*gamma*(c1 + c2*phis)*phi .* dot(nhat, grad(v))" *
              " - theta4*gamma*(c3 + c4*phis)*phis .* dot(nhat, grad(v))" *
              " + boundary(g*v)")

evalInitialConditions()
phi0 = copy(phi.values)
x = g_grid.allnodes[1, :]
y = g_grid.allnodes[2, :]
function exact(t)
    r = sqrt.((x .- xc0 .- ux*t).^2 .+ (y .- yc0 .- uy*t).^2)
    return 0.5 .* (1 .- tanh.((r .- R) ./ (2*eps)))
end
mass0 = sum(g_lump .* phi0[:])
writedlm(joinpath(outdir, "nodes_n$nel.csv"), [x y], ',')

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
    nsteps = ceil(Int, T/(CFL*dx/umax) - 1e-9)
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
    hist = zeros(0, 6)
    status = "ok"; t_end = 0.0; err = NaN; dmass = NaN
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
        err   = maximum(abs.(p .- exact(t_end)))
        dmass = (sum(g_lump .* p) - mass0)/mass0
        hist = vcat(hist, [step t_end err minimum(p) maximum(p) dmass])
        store_frames!(step)
    end
    secs = time() - t0
    p = phi.values[:]
    if status == "ok" && (minimum(p) < -0.01 || maximum(p) > 1.01)
        status = "overshoot"            # finished, but phi left [-0.01, 1.01]
    end
    status == "diverged" && (err = NaN; dmass = NaN)
    writedlm(joinpath(outdir, "hist", name * ".csv"), hist, ',')
    open(joinpath(outdir, "frames", name * ".bin"), "w") do io
        nf = length(ftimes)
        write(io, Int64(g_nnode), Int64(nf), ftimes, frames[:, 1:nf])
    end
    return [name, nel, form, iters, CFL, dt, nsteps, status, t_end, err,
            minimum(p), maximum(p), dmass, round(secs, digits=2)]
end

resfile = joinpath(outdir, "results.csv")
if !isfile(resfile)
    open(resfile, "w") do io
        println(io, "run,nel,form,iterations,CFL,dt,steps,status,t_end,maxerr,phimin,phimax,mass_change,seconds")
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
    @printf("%-28s %-9s t_end %.4f  max error %.4e  phi [%.4f, %.4f]  mass %.1e  %.1f s\n",
            name, row[8], row[9], row[10], row[11], row[12], row[13], row[14])
    flush(stdout)
end

finalizeFinch()
