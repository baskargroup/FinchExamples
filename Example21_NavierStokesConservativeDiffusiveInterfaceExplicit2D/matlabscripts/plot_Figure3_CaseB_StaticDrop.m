% Figure 3 (case B): static drop. Pressure and spurious currents at the end,
% and the pressure along y = 0.5 against the Laplace jump sigma/R = 4.
clc; close all; clear;
setenv('CDI_CASE', 'B');
setenv('CDI_FS', '20');     % readable when the figure fills a slide
S = plot_style();
Q = load_results(S);
sigma = 1; R = 0.25;
A1 = Q.field(Q.A(:,end));
P1 = Q.field(Q.P(:,end));
P1 = P1 - mean(P1(hypot(Q.Xq - 0.5, Q.Yq - 0.5) > 1.6*R));   % far field = 0

fig = S.figure(1500, 480);
tl = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = nexttile(tl); hold(ax, 'on');
S.field(ax, Q.Xq, Q.Yq, P1); colormap(ax, S.cmap); add_colorbar(ax, S, 'Pressure (p)');
contour(ax, Q.Xq, Q.Yq, A1, [0.5 0.5], '-', 'LineColor', 'w', 'LineWidth', 1.5);
style_field_axes(ax, S, Q.Lx, Q.Ly); xticks(ax, 0:0.5:1); yticks(ax, 0:0.5:1);
title(ax, sprintf('Pressure p, t = %.1f', Q.t(end)));

ax = nexttile(tl); hold(ax, 'on');
Sp = Q.field(hypot(Q.UX(:,end), Q.UY(:,end)));
S.field(ax, Q.Xq, Q.Yq, Sp);
colormap(ax, S.cmap); clim(ax, [0 max(Sp(:))]); add_colorbar(ax, S, 'Velocity magnitude |u|');
contour(ax, Q.Xq, Q.Yq, A1, [0.5 0.5], '-', 'LineColor', 'w', 'LineWidth', 1.5);
style_field_axes(ax, S, Q.Lx, Q.Ly); xticks(ax, 0:0.5:1); yticks(ax, 0:0.5:1);
title(ax, sprintf('Velocity |u| (spurious currents), t = %.1f', Q.t(end)));

ax = nexttile(tl); hold(ax, 'on');
row = round(size(Q.Xq, 1)/2);
plot(ax, Q.Xq(row,:), P1(row,:), '-', 'Color', S.blue, 'LineWidth', 2, 'DisplayName', 'Finch');
xs = Q.Xq(row,:); pl = (abs(xs - 0.5) < R) * sigma/R;
plot(ax, xs, pl, '--', 'Color', S.orange, 'LineWidth', 2, 'DisplayName', 'Laplace, \sigma/R');
S.axes(ax); xlabel(ax, 'Position (x), y = 0.5'); ylabel(ax, 'Pressure (p)');
legend(ax, 'Location', 'east', 'Box', 'off', 'TextColor', S.ink);
title(ax, 'Pressure along y = 0.5');
S.save(fig, 'Figure3_CaseB_StaticDrop.pdf');
fprintf('written: %s\n', S.file('Figure3_CaseB_StaticDrop.pdf'));
