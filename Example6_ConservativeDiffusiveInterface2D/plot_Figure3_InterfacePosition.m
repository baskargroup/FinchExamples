% Figure 3: interface position (centroid of phi) vs time
clc; close all; clear;
plot_common;
xc = zeros(nsnap,1);
for k = 1:nsnap
    xc(k) = sum(x .* PHI(:,k)) / sum(PHI(:,k));   % x-centroid of the phase
end
xc_exact = xc0 + ux*times;

fig = S.figure(800, 450);
ax  = axes(fig); hold(ax, 'on');
hN = plot(ax, times, xc, 'o', 'Color', S.blue, 'LineWidth', 1, ...
     'MarkerSize', 6, 'MarkerFaceColor', 'none', 'DisplayName', 'Numerical');
hE = plot(ax, times, xc_exact, '-', 'Color', S.orange, 'LineWidth', 2, ...
     'DisplayName', 'Exact');
S.axes(ax);
xlabel(ax, 'Time'); ylabel(ax, 'Interface position');
xlim(ax, [times(1) times(end)]); ylim(ax, [xc0 xc0 + ux*times(end)]);
grid(ax, 'minor');
ax.MinorGridColor = S.grid; ax.MinorGridAlpha = 0.6; ax.MinorGridLineStyle = ':';
ax.XAxis.MinorTickValues = times(1):0.02:times(end);
ax.YAxis.MinorTickValues = xc0:0.02:(xc0 + ux*times(end));
legend(ax, [hE hN], 'Location', 'northwest', 'Box', 'off', 'TextColor', S.ink2);
S.save(fig, 'Figure3_InterfacePosition.pdf');
fprintf('max interface position error = %.3e\n', max(abs(xc - xc_exact)));
