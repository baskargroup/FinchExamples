% Figure 1: final phase field with numerical and exact interfaces (phi = 0.5)
clc; close all; clear;
plot_common;
T  = times(end);
F.Values = PHI(:,end);
xc = xc0 + ux*T;  yc = yc0 + uy*T;

fig = S.figure(620, 520);
ax  = axes(fig); hold(ax, 'on');
contourf(ax, Xq, Yq, F(Xq, Yq), 20, 'LineColor', 'none');
colormap(ax, S.ramp(64)); clim(ax, [0 1]);
cb = colorbar(ax); cb.Label.String = 'Phase field (\phi)';
cb.Label.Color = S.ink; cb.Color = S.ink2; cb.TickDirection = 'in';
[~, hN] = contour(ax, Xq, Yq, F(Xq, Yq), [0.5 0.5], '-', ...
    'LineColor', S.ink, 'LineWidth', 1.5, ...
    'DisplayName', sprintf('Numerical (t = %g)', T));
hE = plot(ax, xc + R*cos(theta), yc + R*sin(theta), '--', ...
    'Color', S.orange, 'LineWidth', 2, 'DisplayName', sprintf('Exact (t = %g)', T));
style_field_axes(ax, S);
legend(ax, [hE hN], 'Location', 'northwest', 'Box', 'off', 'TextColor', S.ink2);
S.save(fig, 'Figure1_Solution.pdf');

data = readmatrix('solution.csv');
fprintf('max error = %.4e\n', max(abs(data(:,3) - data(:,4))));
