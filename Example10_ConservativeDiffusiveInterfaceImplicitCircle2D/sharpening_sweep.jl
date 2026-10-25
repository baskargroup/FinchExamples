#==============================================================================#
# Sharpening-term study for the semi-implicit 2D CDI scheme of main.jl
# (Finch V2, the scheme in scheme.jl).
#
# The sharpening flux gamma*(phi - phi^2)*nhat is split between the unknown phi
# and a lagged field phis in several ways (scheme.jl):
#   lagged    phis*(1 - phis)              right-hand side only
#   derived   phi*(1 - phis)               as in main.jl
#   reversed  phis*(1 - phi)
#   newton    (1 - 2*phis)*phi + phis^2    Newton about phis
#   extrap    pe*(1 - pe), pe = 2*phio - phioo, one solve per step
#             (BDF2 with an extrapolated explicit term, the usual IMEX choice)
# The first four run with 1, 2 or 3 Picard iterations: iteration 1 lags with
# phio, each later one with the previous iterate. 'lagged' with one iteration
# is the explicit treatment phio*(1 - phio). The normal is always recovered
# from phis, as in main.jl.
#
# Each form runs at six CFL numbers on the mesh given on the command line, with
# eps = dx. dt is rounded down so that T is reached exactly. Measure: the max
# nodal error against the translated tanh circle at T (not an exact solution
# for a curved interface, so the error has a floor set by the mesh). A run
# whose frames contain NaN or leave [-1, 2] is marked diverged (Finch V2 runs
# every step to the end, so the run is not stopped early; t_end is the last
# frame that was still sound).
#
# Writes into output/sharpening_sweep/ (runs already in results.csv are skipped):
#   results.csv          run, nel, form, iterations, CFL, dt, steps, status,
#                        t_end, max error, min phi, max phi, relative mass
#                        change, seconds
#   hist/<run>.csv       step, t, max error, min phi, max phi, mass change, at the frames
#   frames/<run>.bin     phi at 41 instants for the movies: Int64 nnode, Int64 nf,
#                        Float64 times[nf], Float32 phi[nnode, nf]
#   nodes_n<nel>.csv     node coordinates x, y
#
# Usage:  julia sharpening_sweep.jl nel [forms] [iterations] [CFLs]
#   e.g.  julia sharpening_sweep.jl 100 derived 2 0.25
#==============================================================================#

using DelimitedFiles, Printf

include("scheme.jl")

nel = parse(Int, ARGS[1])
listarg(i, f, default) = length(ARGS) >= i ? f.(split(ARGS[i], ",")) : default
FORMS = listarg(2, String, ["extrap", "lagged", "derived", "reversed", "newton"])
ITERS = listarg(3, s -> parse(Int, s), [1, 2, 3])
CFLS  = listarg(4, s -> parse(Float64, s), [2.0, 1.0, 0.5, 0.25, 0.125, 0.0625])

outdir = joinpath(@__DIR__, "output", "sharpening_sweep")
mkpath(joinpath(outdir, "hist")); mkpath(joinpath(outdir, "frames"))

const NF = 40      # movie frames after t = 0
eps_n = 1.0/nel
function exact(x, y, t)
    r = sqrt.((x .- XC0 .- UX*t).^2 .+ (y .- YC0 .- UY*t).^2)
    return 0.5 .* (1 .- tanh.((r .- R) ./ (2*eps_n)))
end

# One run to t = T. Returns the results row; writes the history and frames.
function sweep_case(form, iters, CFL, name)
    t0 = time()
    r = run_case(CFL, iters; form = form, nel = nel, snapshots = NF)
    secs = time() - t0
    x, y = r.x, r.y
    isfile(joinpath(outdir, "nodes_n$nel.csv")) || writedlm(joinpath(outdir, "nodes_n$nel.csv"), [x y], ',')
    # Frames up to the first unsound one
    sound(p) = all(isfinite, p) && maximum(p) <= 2 && minimum(p) >= -1
    nf = r.ok ? something(findfirst(p -> !sound(p), r.snapshots), length(r.snapshots) + 1) - 1 : 0
    status = (r.ok && nf == length(r.snapshots)) ? "ok" : "diverged"
    frames = r.snapshots[1:nf]; ftimes = r.snap_times[1:nf]
    mass0 = r.ok ? r.mass[1] : NaN
    hist = zeros(0, 6)
    for (k, p) in enumerate(frames)
        step = round(Int, ftimes[k]/r.dt)
        hist = vcat(hist, [step ftimes[k] maximum(abs.(p .- exact(x, y, ftimes[k]))) minimum(p) maximum(p) (r.mass[step + 1] - mass0)/mass0])
    end
    t_end = nf > 0 ? ftimes[end] : 0.0
    p = status == "ok" ? r.phi : (nf > 0 ? frames[end] : fill(NaN, length(x)))
    err   = status == "ok" ? maximum(abs.(p .- exact(x, y, T))) : NaN
    dmass = status == "ok" ? (r.mass[end] - mass0)/mass0 : NaN
    if status == "ok" && (minimum(p) < -0.01 || maximum(p) > 1.01)
        status = "overshoot"            # finished, but phi left [-0.01, 1.01]
    end
    writedlm(joinpath(outdir, "hist", name * ".csv"), hist, ',')
    open(joinpath(outdir, "frames", name * ".bin"), "w") do io
        write(io, Int64(length(x)), Int64(nf), ftimes,
              nf > 0 ? Float32.(reduce(hcat, frames)) : zeros(Float32, 0))
    end
    return [name, nel, form, iters, CFL, r.dt, r.nsteps, status, t_end, err,
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
    row = sweep_case(form, iters, CFL, name)
    open(resfile, "a") do io
        println(io, join(row, ","))
    end
    @printf("%-28s %-9s t_end %.4f  max error %.4e  phi [%.4f, %.4f]  mass %.1e  %.1f s\n",
            name, row[8], row[9], row[10], row[11], row[12], row[13], row[14])
    flush(stdout)
end
