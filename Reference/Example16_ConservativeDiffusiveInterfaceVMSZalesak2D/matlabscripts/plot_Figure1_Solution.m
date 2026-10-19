% Figure 1: phase field after one revolution, with the phi = 0.5 contour of the
% final field and of the initial field (= exact solution)
clc; close all; clear;
plot_common;
T  = times(end);
F.Values = PHI(:,end);
Zf = F(Xq, Yq);
F.Values = PHI(:,1);
Z0 = F(Xq, Yq);

fig = S.figure(620, 520);
ax  = axes(fig); hold(ax, 'on');
contourf(ax, Xq, Yq, Zf, 20, 'LineColor', 'none');
colormap(ax, S.ramp(64)); clim(ax, [0 1]);
cb = colorbar(ax); cb.Label.String = 'Phase field (\phi)';
cb.Label.Color = S.ink; cb.Color = S.ink2; cb.TickDirection = 'in';
[~, hN] = contour(ax, Xq, Yq, Zf, [0.5 0.5], '-', ...
    'LineColor', S.ink, 'LineWidth', 1.5, ...
    'DisplayName', sprintf('Numerical (t = %.3f)', T));
[~, hE] = contour(ax, Xq, Yq, Z0, [0.5 0.5], '--', ...
    'LineColor', S.orange, 'LineWidth', 2, 'DisplayName', 'Exact (initial field)');
style_field_axes(ax, S);
legend(ax, [hE hN], 'Location', 'southwest', 'Box', 'off', 'TextColor', S.ink2);
S.save(fig, 'Figure1_Solution.pdf');

data = S.load('solution.csv');
fprintf('max error = %.4e\n', max(abs(data(:,3) - data(:,4))));
