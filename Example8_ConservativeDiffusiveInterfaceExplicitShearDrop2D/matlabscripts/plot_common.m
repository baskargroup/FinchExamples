% Shared loading for the Example 8 plot scripts (called by each script).
% solution.csv:            x, y, phi (t = T), phi (initial), phi (t = T/2)
% solution_timeseries.csv: columns 1-2 = x, y; columns 3.. = phi at snapshot k
% snapshot_times.csv:      time of each snapshot
S     = plot_style();
ts    = S.load('solution_timeseries.csv');
times = S.load('snapshot_times.csv');
x     = ts(:,1);
y     = ts(:,2);
PHI   = ts(:,3:end);
nsnap = size(PHI,2);

% Problem parameters (must match main.jl)
xc  = 0.5;  yc = 0.75;       % drop centre
R   = 0.15;                  % drop radius
Tp  = 4;                     % period of the velocity field
dx  = min(diff(unique(x)));  % mesh size

% Regular grid for contour plotting, at least twice as fine as the mesh
nx = max(400, 2*round(1/dx)); ny = nx;
xq = linspace(min(x), max(x), nx);
yq = linspace(min(y), max(y), ny);
[Xq, Yq] = meshgrid(xq, yq);
F = scatteredInterpolant(x, y, PHI(:,1), 'linear', 'none');

% Outline of the initial drop (= exact solution at t = T), for overlays
th = linspace(0, 2*pi, 400);
outline = [xc + R*cos(th);  yc + R*sin(th)];
