% Figure 3: interface position (phi = 0.5 crossing) vs time
clc; close all; clear;
% solution_timeseries.csv: column 1 = x, columns 2..N+1 = phi at snapshot k
% snapshot_times.csv:      time of each snapshot
S     = plot_style();
ts    = S.load('solution_timeseries.csv');
times = S.load('snapshot_times.csv');
x     = ts(:,1);
PHI   = ts(:,2:end);
nsnap = size(PHI,2);

% Problem parameters (must match main.jl)
u   = 1.0;
x0  = 0.25;
eps = 1/200;

% Interface position = phi = 0.5 crossing, linearly interpolated between the
% last node below 0.5 and the first node above it (sub-element accuracy).
xi = zeros(nsnap,1);
for k = 1:nsnap
    j = find(PHI(:,k) >= 0.5, 1);            % first node with phi >= 0.5
    p0 = PHI(j-1,k); p1 = PHI(j,k);
    xi(k) = x(j-1) + (0.5 - p0)/(p1 - p0) * (x(j) - x(j-1));
end
xi_exact = x0 + u*times;

fig3 = S.figure(800, 450);
ax = axes(fig3); hold(ax, 'on');
hN = plot(ax, times, xi, 'o', 'Color', S.blue, 'LineWidth', 1, ...
     'MarkerSize', 6, 'MarkerFaceColor', 'none', 'DisplayName', 'Numerical');
hE = plot(ax, times, xi_exact, '-', 'Color', S.orange, 'LineWidth', 2, ...
     'DisplayName', 'Exact');
S.axes(ax);
xlabel(ax, 'Time'); ylabel(ax, 'Interface position');
xlim(ax, [times(1) times(end)]); ylim(ax, [0.25 0.5]);
grid(ax, 'minor');                       % major + minor grid
ax.MinorGridColor = S.grid; ax.MinorGridAlpha = 0.6; ax.MinorGridLineStyle = ':';
ax.XAxis.MinorTickValues = 0:0.0125:0.25;
ax.YAxis.MinorTickValues = 0.25:0.0125:0.5;
legend(ax, [hE hN], 'Location', 'northwest', 'Box', 'off', 'TextColor', S.ink2);
S.save(fig3, 'Figure3_InterfacePosition.pdf');
fprintf('max interface position error = %.3e\n', max(abs(xi - xi_exact)));
