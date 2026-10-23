% Shared loading for the Example 13 plot scripts (called by each script).
% 1D files: 1d_fd.csv      x, phi, forward, backward, central, exact
%           1d_nodal.csv   x, phi, nodal gradient, exact
%           1d_elem.csv    element, node, x, element gradient (one row per
%                          element and local node)
% 2D files: 2d_fd.csv      x, y, phi, forward (x,y), backward (x,y), central (x,y)
%           2d_nodal.csv   x, y, phi, nodal gradient (x,y), exact normal (x,y)
%           2d_elem.csv    element, node, x, y, element gradient (x,y)
S = plot_style();

% Problem parameters (must match main.jl)
h   = 1/40;                  % node spacing = element size
eps = h;                     % interface thickness
x0  = 0.5 + 0.3*h;           % 1D interface position
xc  = 0.5 + 0.3*h;           % 2D disk centre and radius
yc  = 0.5 + 0.2*h;
R   = 0.25;

% Colours (validated categorical slots 1-3, fixed order)
C.central  = S.blue;                     % centred: FD central, FEM nodal
C.exact    = S.orange;
C.onesided = [ 27 175 122]/255;          % one-sided: FD forward/backward, FEM element-wise

% Angle (degrees) from the exact normal, positive counter-clockwise; the
% exact normal points to the disk centre
angerr = @(gx, gy, x, y) atan2d(gx.*(-(y-yc)) - gy.*(-(x-xc)), ...
                                gx.*(-(x-xc)) + gy.*(-(y-yc)));
polar  = @(x, y) mod(atan2d(y - yc, x - xc), 360);
inband = @(p) p > 0.05 & p < 0.95;       % the interface band used for errors
