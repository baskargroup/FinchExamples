% Figure 1 (case A): drop advection with rho1/rho2 = 10. Consistent momentum
% transport (theta23 = 1) against the inconsistent one (theta23 = 0): the
% drop after one period with the velocity error, and the error in time.
% Needs the runs output_caseA and output_caseA_theta23_0; the inconsistent
% run blows up at t = 0.94, so it is run to t = 0.9 (see README).
clc; close all; clear;
setenv('CDI_CASE', 'A');
setenv('CDI_FS', '20');     % readable when the figure fills a slide
S = plot_style();
runs   = {'output_caseA', 'output_caseA_theta23_0'};
labels = {'With consistent term', 'Without consistent term'};
R = cellfun(@(f) load_results(S, f), runs, 'UniformOutput', false);


fig = S.figure(1000, 1150);
tl = tiledlayout(fig, 3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for r = 1:2
    Q = R{r};
    A0 = Q.field(Q.A(:,1)); A1 = Q.field(Q.A(:,end));
    ax = nexttile(tl); hold(ax, 'on');
    S.field(ax, Q.Xq, Q.Yq, A1);
    colormap(ax, S.cmap); clim(ax, [0 1]); add_colorbar(ax, S, 'Phase field (\phi)');
    contour(ax, Q.Xq, Q.Yq, A0, [0.5 0.5], '--', 'LineColor', 'w', 'LineWidth', 1.5);
    contour(ax, Q.Xq, Q.Yq, A1, [0.5 0.5], '-', 'LineColor', S.ink, 'LineWidth', 1.5);
    style_field_axes(ax, S, Q.Lx, Q.Ly); xticks(ax, 0:0.5:1); yticks(ax, 0:0.5:1);
    title(ax, sprintf('%s, t = %.2f', labels{r}, Q.t(end)));

    ax = nexttile(tl); hold(ax, 'on');
    E = Q.field(hypot(Q.UX(:,end) - 1, Q.UY(:,end)));
    S.field(ax, Q.Xq, Q.Yq, E);
    colormap(ax, S.cmap); clim(ax, [0 max(max(E(:)), eps)]);   % own range: the errors differ by 10^4
    add_colorbar(ax, S, '|u - (1, 0)|');
    contour(ax, Q.Xq, Q.Yq, A1, [0.5 0.5], '-', 'LineColor', 'w', 'LineWidth', 1.5);
    style_field_axes(ax, S, Q.Lx, Q.Ly); xticks(ax, 0:0.5:1); yticks(ax, 0:0.5:1);
    title(ax, 'Velocity error');
end

ax = nexttile(tl, [1 2]); hold(ax, 'on');
cols = {S.blue, S.orange};
for r = 1:2
    h = R{r}.h;
    semilogy(ax, h(:,1), max(h(:,10), 1e-16), '-', 'Color', cols{r}, 'LineWidth', 2, ...
             'DisplayName', labels{r});
end
set(ax, 'YScale', 'log');
S.axes(ax);
xlabel(ax, 'Time (t)'); ylabel(ax, '||u - (1, 0)||_2');
legend(ax, 'Location', 'southeast', 'Box', 'off', 'TextColor', S.ink);
title(ax, 'Velocity error in time (exact solution: u = (1, 0) for all t)');

S.save(fig, 'Figure1_CaseA_DropAdvection.pdf');
fprintf('written: %s\n', S.file('Figure1_CaseA_DropAdvection.pdf'));
