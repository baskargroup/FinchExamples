% Figure 8 (case D): dam break, Martin and Moyce (1952) column of width a
% and height 2a. Left: the water (phi = 0.5 contour on the phase field) at
% four times. Right: surge front position Z = x/a and column height
% H = h/(2a) against T = t sqrt(2g/a), with the experimental front from
% reference_martin_moyce/front_n2.txt.
clc; close all; clear;
setenv('CDI_CASE', 'D');
setenv('CDI_FS', '18');
S = plot_style();
Q = load_results(S);
a = 0.25; g = 1;
F = readmatrix(S.file('front.csv'));          % t, x_front, h_left, T, Z, H
E = readmatrix(fullfile(fileparts(mfilename('fullpath')), '..', 'reference_martin_moyce', 'front_n2.txt'), 'FileType', 'text', 'CommentStyle', '#');

fig = S.figure(1500, 760);
tl = tiledlayout(fig, 2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
Tshow = [0 1 2 3];
for k = 1:4
    [~, j] = min(abs(Q.t * sqrt(2 * g / a) - Tshow(k)));
    ax = nexttile(tl, k); hold(ax, 'on');
    A = Q.field(Q.A(:,j));
    S.field(ax, Q.Xq, Q.Yq, A); colormap(ax, S.cmap); clim(ax, [0 1]);
    contour(ax, Q.Xq, Q.Yq, A, [0.5 0.5], '-', 'LineColor', 'w', 'LineWidth', 1.5);
    S.axes(ax); axis(ax, 'equal'); xlim(ax, [0 Q.Lx]); ylim(ax, [0 Q.Ly]); ax.Layer = 'top';
    xticks(ax, 0:0.5:1); yticks(ax, 0:0.25:0.75); grid(ax, 'off');
    if k == 1, ylabel(ax, 'Position (y)'); end
    xlabel(ax, 'Position (x)');
    title(ax, sprintf('T = %.1f', Q.t(j) * sqrt(2 * g / a)));
end
ax = nexttile(tl, 5, [1 2]); hold(ax, 'on');
plot(ax, F(:,4), F(:,5), '-', 'Color', S.blue, 'LineWidth', 2.5, 'DisplayName', 'Finch, 128 x 96');
plot(ax, E(:,1), E(:,2), 'o', 'Color', S.orange, 'MarkerFaceColor', S.orange, 'MarkerSize', 7, ...
     'DisplayName', 'Martin and Moyce (1952)');
S.axes(ax); xlim(ax, [0 3.5]); ylim(ax, [1 4]);
xlabel(ax, 'T = t (2g/a)^{1/2}'); ylabel(ax, 'Front position x / a');
legend(ax, 'Location', 'northwest', 'Box', 'off', 'TextColor', S.ink);
ax = nexttile(tl, 7, [1 2]); hold(ax, 'on');
plot(ax, F(:,4), F(:,6), '-', 'Color', S.blue, 'LineWidth', 2.5, 'DisplayName', 'Finch, 128 x 96');
S.axes(ax); xlim(ax, [0 max(F(:,4))]); ylim(ax, [0 1.05]);
xlabel(ax, 'T = t (2g/a)^{1/2}'); ylabel(ax, 'Column height h / 2a');
S.save(fig, 'Figure8_CaseD_DamBreak.pdf');
fprintf('written: %s\n', S.file('Figure8_CaseD_DamBreak.pdf'));
