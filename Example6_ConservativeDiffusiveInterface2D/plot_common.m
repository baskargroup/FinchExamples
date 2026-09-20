% Shared loading for the Example 6 plot scripts (called by each script).
% solution.csv:            x, y, phi (numerical), phi (exact)  [final time]
% solution_timeseries.csv: columns 1-2 = x, y; columns 3.. = phi at snapshot k
% snapshot_times.csv:      time of each snapshot
S     = plot_style();
ts    = readmatrix('solution_timeseries.csv');
times = readmatrix('snapshot_times.csv');
x     = ts(:,1);
y     = ts(:,2);
PHI   = ts(:,3:end);
nsnap = size(PHI,2);

% Problem parameters (must match main.jl)
ux  = 1.0;  uy = 0.0;        % advection velocity
xc0 = 0.3;  yc0 = 0.5;       % initial circle center
R   = 0.15;                  % circle radius
dx  = min(diff(unique(x)));  % mesh size (eps = dx in main.jl)
eps = dx;

% Regular grid for contour plotting
nx = 300; ny = 300;
xq = linspace(min(x), max(x), nx);
yq = linspace(min(y), max(y), ny);
[Xq, Yq] = meshgrid(xq, yq);
F = scatteredInterpolant(x, y, PHI(:,1), 'linear', 'none');
theta = linspace(0, 2*pi, 400);
