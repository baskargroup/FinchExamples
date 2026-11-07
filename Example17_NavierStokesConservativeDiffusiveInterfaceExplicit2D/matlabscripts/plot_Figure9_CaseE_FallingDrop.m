% Figure 9 (case E): a water drop falling in air into a pool. Top: the phase
% field at six times. Bottom: centroid height and mean vertical velocity of
% the water above the pool against free fall, y = y0 - g t^2/2, v = -g t,
% until the drop touches the pool (drop.csv from main.jl).
clc; close all; clear;
setenv('CDI_CASE', 'E');
setenv('CDI_FS', '18');
S = plot_style();
Q = load_results(S);
g = 1; y0 = 1.5; R = 0.1; pool = 0.5;
% t, y_c, v_c of the water above the pool (y > pool + 4 eps), from the snapshots
eps = 1 / 96; above = Q.y > pool + 4 * eps;
w = Q.A(above, :); D = [Q.t(:), (w' * Q.y(above)) ./ max(sum(w, 1)', 1e-300), sum(w .* Q.UY(above, :), 1)' ./ max(sum(w, 1)', 1e-300)];
t_imp = sqrt(2 * (y0 - R - pool) / g);        % free-fall impact time
m = D(:,1) <= t_imp;

fig = S.figure(1500, 900);
tl = tiledlayout(fig, 2, 6, 'TileSpacing', 'compact', 'Padding', 'compact');
tshow = [0 0.7 1.33 1.5 1.8 2.4];
for k = 1:6
    [~, j] = min(abs(Q.t - tshow(k)));
    ax = nexttile(tl, k); hold(ax, 'on');
    A = Q.field(Q.A(:,j));
    S.field(ax, Q.Xq, Q.Yq, A); colormap(ax, S.cmap); clim(ax, [0 1]);
    contour(ax, Q.Xq, Q.Yq, A, [0.5 0.5], '-', 'LineColor', 'w', 'LineWidth', 1.5);
    S.axes(ax); axis(ax, 'equal'); xlim(ax, [0 Q.Lx]); ylim(ax, [0 Q.Ly]); ax.Layer = 'top'; grid(ax, 'off');
    xticks(ax, 0:0.5:1); yticks(ax, 0:0.5:2);
    if k == 1, ylabel(ax, 'Position (y)'); else, yticklabels(ax, []); end
    xlabel(ax, 'Position (x)');
    title(ax, sprintf('t = %.2f', Q.t(j)));
end
ax = nexttile(tl, 7, [1 3]); hold(ax, 'on');
plot(ax, D(m,1), D(m,2), '-', 'Color', S.blue, 'LineWidth', 2.5, 'DisplayName', 'Finch, 96 x 192');
plot(ax, D(m,1), y0 - g * D(m,1).^2 / 2, '--', 'Color', S.orange, 'LineWidth', 2, 'DisplayName', 'free fall, y_0 - g t^2/2');
S.axes(ax); xlim(ax, [0 t_imp]); xlabel(ax, 'Time (t)'); ylabel(ax, 'Drop centroid (y)');
legend(ax, 'Location', 'southwest', 'Box', 'off', 'TextColor', S.ink);
ax = nexttile(tl, 10, [1 3]); hold(ax, 'on');
plot(ax, D(m,1), -D(m,3), '-', 'Color', S.blue, 'LineWidth', 2.5, 'DisplayName', 'Finch, 96 x 192');
plot(ax, D(m,1), g * D(m,1), '--', 'Color', S.orange, 'LineWidth', 2, 'DisplayName', 'free fall, g t');
S.axes(ax); xlim(ax, [0 t_imp]); xlabel(ax, 'Time (t)'); ylabel(ax, 'Drop speed');
legend(ax, 'Location', 'northwest', 'Box', 'off', 'TextColor', S.ink);
S.save(fig, 'Figure9_CaseE_FallingDrop.pdf');
fprintf('written: %s\n', S.file('Figure9_CaseE_FallingDrop.pdf'));
