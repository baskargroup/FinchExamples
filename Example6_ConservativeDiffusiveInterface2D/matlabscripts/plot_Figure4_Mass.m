% Figure 4: total mass (integral of phi) relative to its initial value vs time
clc; close all; clear;
plot_common;
% Lumped finite element weights: the lumped mass entry of a node is the integral
% of its basis function, which the boundary truncates to half at an edge node and
% a quarter at a corner. This keeps the sum exact once phi reaches the wall.
w = ones(size(x));
w(abs(x) < 1e-12 | abs(x-1) < 1e-12) = 0.5 * w(abs(x) < 1e-12 | abs(x-1) < 1e-12);
w(abs(y) < 1e-12 | abs(y-1) < 1e-12) = 0.5 * w(abs(y) < 1e-12 | abs(y-1) < 1e-12);
mass = (w' * PHI)' * dx * dx;
rel  = mass / mass(1) - 1;

fig = S.figure(800, 450);
ax  = axes(fig); hold(ax, 'on');
yline(ax, 0, '-', 'Color', S.orange, 'LineWidth', 2);
plot(ax, times, rel, 'o-', 'Color', S.blue, 'LineWidth', 1.2, ...
     'MarkerSize', 6, 'MarkerFaceColor', 'none');
S.axes(ax);
xlabel(ax, 'Time'); ylabel(ax, 'Relative mass change');
xlim(ax, [times(1) times(end)]);
m = max(abs(rel)); if m == 0, m = 1e-12; end
ylim(ax, 1.3*[-m m]);
grid(ax, 'minor');
ax.MinorGridColor = S.grid; ax.MinorGridAlpha = 0.6; ax.MinorGridLineStyle = ':';
ax.XAxis.MinorTickValues = times(1):0.02:times(end);
S.save(fig, 'Figure4_Mass.pdf');
fprintf('relative mass change at final time = %.3e\n', rel(end));
