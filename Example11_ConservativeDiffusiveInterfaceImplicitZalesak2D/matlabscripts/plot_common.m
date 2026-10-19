% Shared loading for the Example 11 plot scripts (called by each script).
% solution.csv:            x, y, phi (after one revolution), phi (initial)
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
xc  = 0.5;  yc = 0.75;       % disk centre
R   = 0.15;                  % disk radius
sw  = 0.05;                  % slot width
sy1 = yc - R;  sy2 = yc + 0.10;   % slot mouth (disk bottom) and slot top; main.jl's
                                  % distance-function rectangle extends below sy1
dx  = min(diff(unique(x)));  % mesh size (eps = dx in main.jl)
eps = dx;

% Regular grid for contour plotting
nx = 400; ny = 400;
xq = linspace(min(x), max(x), nx);
yq = linspace(min(y), max(y), ny);
[Xq, Yq] = meshgrid(xq, yq);
F = scatteredInterpolant(x, y, PHI(:,1), 'linear', 'none');

% Outline of the exact slotted disk (for overlays): disk arc from the right
% slot edge round to the left slot edge, then the slot.
a1 = asin((sw/2)/R);
th = linspace(-pi/2 + a1, 3*pi/2 - a1, 400);
outline = [xc + R*cos(th);  yc + R*sin(th)];
outline = [outline, [xc - sw/2, xc - sw/2, xc + sw/2, xc + sw/2; ...
                     yc - R*cos(a1), sy2, sy2, yc - R*cos(a1)]];
outline = [outline, outline(:,1)];
