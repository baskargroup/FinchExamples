#==============================================================================#
# Interface normal of the CDI model: finite differences vs finite elements.
#
# The CDI sharpening term needs the unit normal nhat = grad(phi)/|grad(phi)|
# of a tanh profile that is only a few cells wide (eps = h). This example
# computes that gradient and normal on a fixed field, with no time stepping,
# by five methods on the SAME nodes (spacing h):
#
#   finite differences on the nodal values
#     FD forward     (phi_{i+1} - phi_i)/h
#     FD backward    (phi_i - phi_{i-1})/h
#     FD central     (phi_{i+1} - phi_{i-1})/(2h)   the normal of Mirjalili et
#                    al., JCP 401 (2020), Eq. 11; each component from its own
#                    direction only
#   finite elements, linear Lagrange basis (as in Examples 5-12)
#     element-wise   gradient of the element's own interpolant, evaluated at
#                    the element's nodes: a node shared by several elements
#                    gets several different values (the jump)
#     nodal          lumped L2 projection of the element gradients onto the
#                    nodal basis, then normalized at the nodes (recover_nhat!
#                    of Examples 5-12): one continuous value per node
#
# The analogy being demonstrated, on a uniform grid:
#   - linear element-wise gradient at a node = the one-sided difference of
#     the element it is taken from: the element to the right (1D) or to the
#     upper right (2D) gives FD forward, to the left / lower left gives FD
#     backward. The element-wise normal is a one-sided difference, and it
#     jumps between elements by forward - backward.
#   - linear nodal gradient = FD central exactly in 1D; in 2D it is the
#     central difference averaged over the neighbouring rows with weights
#     1/6, 2/3, 1/6 (a 9-point stencil, centred like the paper's 3-point one;
#     the weights come from integrating the bilinear element gradient, which
#     varies linearly across the element, against the basis).
# Both identities are checked numerically below (printed differences ~1e-14).
#
# Higher-order elements do not remove the jump, only shrink it: with quadratic
# elements on the same nodes (an earlier version of this example) the angle
# error on the disk was 8.9 deg element-wise and 4.9 deg nodal (worst at the
# edge-midpoint nodes, where the lumped projection averages unevenly), against
# 20 deg linear element-wise and 0.07 deg linear nodal.
#
# Test fields (exact gradient known):
#   1D: phi = 0.5*(1 - tanh((x - x0)/(2 eps))) on [0, 1]
#   2D: a disk, phi = 0.5*(1 - tanh((r - R)/(2 eps))); the exact normal points
#       to the centre (phi = 1 inside), so the angle error of each method is
#       measured against it.
# The interface is placed off the nodes (x0 and the centre are shifted by a
# fraction of h) so that no method is helped by symmetry.
#
# Run as (one dimension per Finch session; run.sh runs both):
#   CDI_DIM=1|2  julia --project=../Finch_patched main.jl
#==============================================================================#

using Finch
using DelimitedFiles

const DIM   = parse(Int, get(ENV, "CDI_DIM", "2"))
tag = string(DIM) * "d"

initFinch("normal" * tag);
useLog("normal" * tag * "log", dir="output", level=1)

domain(DIM)
functionSpace(order=1)

const NNODE = 40                  # elements (node intervals) per direction
const h     = 1.0/NNODE
mesh(DIM == 1 ? LINEMESH : QUADMESH, elsperdim=NNODE)

phi = variable("phi")             # Finch needs one variable to build the problem

# Test field
const eps = h                     # interface thickness, as in the examples
const x0  = 0.5 + 0.3*h           # 1D interface position, off the nodes
const xc  = 0.5 + 0.3*h           # 2D disk centre, off the nodes
const yc  = 0.5 + 0.2*h
const R   = 0.25

phi1d(x)    = 0.5*(1 - tanh((x - x0)/(2*eps)))
dphi1d(x)   = -0.25/eps * sech((x - x0)/(2*eps))^2
rad(x, y)   = sqrt((x - xc)^2 + (y - yc)^2)
phi2d(x, y) = 0.5*(1 - tanh((rad(x, y) - R)/(2*eps)))
# exact unit normal grad(phi)/|grad(phi)|: towards the centre
nexact(x, y) = (-(x - xc)/rad(x, y), -(y - yc)/rad(x, y))

const grid  = Finch.finch_state.grid_data
const refel = Finch.finch_state.refel
const geo   = Finch.finch_state.geo_factors
const Np    = refel.Np
const Nqp   = refel.Nqp
const nodes = grid.allnodes
const nn    = size(nodes, 2)
const ne    = size(grid.loc2glb, 2)

# Nodal values of the test field
p = DIM == 1 ? [phi1d(nodes[1, k]) for k = 1:nn] :
               [phi2d(nodes[1, k], nodes[2, k]) for k = 1:nn]

# Jacobian entries of element e (constant per element on this mesh); in 1D
# Finch leaves the y-parts empty
function jac(e)
    J = geo.J[e]
    DIM == 1 && return J.rx[1], 0.0, 0.0, 0.0
    return J.rx[1], J.ry[1], J.sx[1], J.sy[1]
end

# Element-wise gradient at each element's own nodes: Ddr, Dds are the
# derivative matrices at the element nodes; the Jacobian is constant per
# element on this mesh. One row per (element, local node).
function elementwise_gradient(p)
    rows = zeros(ne*Np, 2 + 2*DIM)          # e, node, x(,y), gx(,gy)
    k = 0
    for e = 1:ne
        rx, ry, sx, sy = jac(e)
        pe = p[grid.loc2glb[:, e]]
        dr = refel.Ddr * pe
        ds = DIM == 2 ? refel.Dds * pe : zeros(Np)
        for i = 1:Np
            k += 1
            g = grid.loc2glb[i, e]
            if DIM == 1
                rows[k, :] = [e, g, nodes[1, g], dr[i]*rx]
            else
                gx = dr[i]*rx + ds[i]*sx
                gy = dr[i]*ry + ds[i]*sy
                rows[k, :] = [e, g, nodes[1, g], nodes[2, g], gx, gy]
            end
        end
    end
    return rows
end

# Nodal gradient: lumped L2 projection of the element gradient onto the
# nodal basis (the recover_nhat! of Examples 5-12, in 1D or 2D):
#     g_l = sum_e int_e (dphi_h/dx) psi_l / sum_e int_e psi_l
function nodal_gradient(p)
    lump = zeros(nn)
    gx = zeros(nn); gy = zeros(nn)
    for e = 1:ne
        rx, ry, sx, sy = jac(e)
        detj = geo.detJ[1, e]
        pe   = p[grid.loc2glb[:, e]]
        dr = refel.Qr * pe                   # d/dr at the quadrature points
        ds = DIM == 2 ? refel.Qs * pe : zeros(Nqp)
        for i = 1:Np
            l = grid.loc2glb[i, e]
            for q = 1:Nqp
                w = refel.wg[q] * refel.Q[q, i] * detj
                lump[l] += w
                gx[l]   += w * (dr[q]*rx + ds[q]*sx)
                gy[l]   += w * (dr[q]*ry + ds[q]*sy)
            end
        end
    end
    return gx ./ lump, gy ./ lump
end

ew = elementwise_gradient(p)
ngx, ngy = nodal_gradient(p)

# Finite differences on the nodal grid (interior nodes; boundary rows NaN).
# Nodes are sorted onto an (NNODE+1)^DIM index grid first, since Finch's
# node numbering is not lexicographic.
idx(c) = round(Int, c/h) + 1
if DIM == 1
    P = fill(NaN, NNODE + 1)
    for k = 1:nn; P[idx(nodes[1, k])] = p[k]; end
    fwd = fill(NaN, nn); bwd = fill(NaN, nn); cen = fill(NaN, nn)
    for k = 1:nn
        i = idx(nodes[1, k])
        (1 < i < NNODE + 1) || continue
        fwd[k] = (P[i+1] - P[i])/h
        bwd[k] = (P[i] - P[i-1])/h
        cen[k] = (P[i+1] - P[i-1])/(2h)
    end
else
    P = fill(NaN, NNODE + 1, NNODE + 1)
    for k = 1:nn; P[idx(nodes[1, k]), idx(nodes[2, k])] = p[k]; end
    fwd = fill(NaN, nn, 2); bwd = fill(NaN, nn, 2); cen = fill(NaN, nn, 2)
    for k = 1:nn
        i = idx(nodes[1, k]); j = idx(nodes[2, k])
        (1 < i < NNODE + 1 && 1 < j < NNODE + 1) || continue
        fwd[k, :] = [(P[i+1, j] - P[i, j])/h,     (P[i, j+1] - P[i, j])/h]
        bwd[k, :] = [(P[i, j] - P[i-1, j])/h,     (P[i, j] - P[i, j-1])/h]
        cen[k, :] = [(P[i+1, j] - P[i-1, j])/(2h), (P[i, j+1] - P[i, j-1])/(2h)]
    end
end

# Identity checks: nodal = central (1D), element-wise =
# forward/backward of the element to the right/upper right and left/lower left
function check_identities()
    if DIM == 1
        ok = .!isnan.(cen)
        println("max |P1 nodal - FD central|  = ", maximum(abs.(ngx[ok] .- cen[ok])))
        dfwd = 0.0; dbwd = 0.0
        for r = 1:size(ew, 1)
            k = Int(ew[r, 2]); isnan(fwd[k]) && continue
            e = Int(ew[r, 1]); xl = minimum(nodes[1, grid.loc2glb[:, e]])
            if nodes[1, k] ≈ xl          # node is the element's left end
                dfwd = max(dfwd, abs(ew[r, 4] - fwd[k]))
            else
                dbwd = max(dbwd, abs(ew[r, 4] - bwd[k]))
            end
        end
        println("max |P1 element-wise - FD forward|  (element to the right) = ", dfwd)
        println("max |P1 element-wise - FD backward| (element to the left)  = ", dbwd)
    else
        dfwd = 0.0; dbwd = 0.0
        for r = 1:size(ew, 1)
            k = Int(ew[r, 2]); isnan(fwd[k, 1]) && continue
            e = Int(ew[r, 1])
            xl = minimum(nodes[1, grid.loc2glb[:, e]]); yl = minimum(nodes[2, grid.loc2glb[:, e]])
            if nodes[1, k] ≈ xl && nodes[2, k] ≈ yl            # lower-left corner
                dfwd = max(dfwd, maximum(abs.(ew[r, 5:6] .- fwd[k, :])))
            elseif !(nodes[1, k] ≈ xl) && !(nodes[2, k] ≈ yl)  # upper-right corner
                dbwd = max(dbwd, maximum(abs.(ew[r, 5:6] .- bwd[k, :])))
            end
        end
        println("max |Q1 element-wise - FD forward|  (element to the upper right) = ", dfwd)
        println("max |Q1 element-wise - FD backward| (element to the lower left)  = ", dbwd)
        # nodal Q1 = central difference averaged over neighbouring rows (1/6, 2/3, 1/6)
        dsm = 0.0
        for k = 1:nn
            i = idx(nodes[1, k]); j = idx(nodes[2, k])
            (2 < i < NNODE && 2 < j < NNODE) || continue
            sx = ((P[i+1, j-1] - P[i-1, j-1])/6 + 2*(P[i+1, j] - P[i-1, j])/3 +
                  (P[i+1, j+1] - P[i-1, j+1])/6)/(2h)
            sy = ((P[i-1, j+1] - P[i-1, j-1])/6 + 2*(P[i, j+1] - P[i, j-1])/3 +
                  (P[i+1, j+1] - P[i+1, j-1])/6)/(2h)
            dsm = max(dsm, abs(ngx[k] - sx), abs(ngy[k] - sy))
        end
        println("max |Q1 nodal - row-averaged FD central| = ", dsm)
    end
end
check_identities()

# Output
if DIM == 1
    x = nodes[1, :]
    writedlm("output/1d_nodal.csv", [x p ngx dphi1d.(x)], ',')
    writedlm("output/1d_elem.csv", ew, ',')
    writedlm("output/1d_fd.csv", [x p fwd bwd cen dphi1d.(x)], ',')
else
    x = nodes[1, :]; y = nodes[2, :]
    ex = [nexact(x[k], y[k]) for k = 1:nn]
    writedlm("output/2d_nodal.csv", [x y p ngx ngy first.(ex) last.(ex)], ',')
    writedlm("output/2d_elem.csv", ew, ',')
    writedlm("output/2d_fd.csv",
             [x y p fwd[:, 1] fwd[:, 2] bwd[:, 1] bwd[:, 2] cen[:, 1] cen[:, 2]], ',')
end

finalizeFinch()
println("wrote output for $(DIM)D: $(nn) nodes, $(ne) elements")
