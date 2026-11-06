% Figure 2 (case C): rising bubble, Hysing et al. (2009) test case 1.
% Left: the bubble at t = 0 and t = 3 with the velocity. Right: centroid,
% rise velocity and circularity in time, with the benchmark values
% (max rise velocity 0.2417 at t = 0.9213, min circularity 0.9013 at
% t = 1.9041, centroid 1.0813 at t = 3).
clc; close all; clear;
setenv('CDI_CASE', 'C');
setenv('CDI_FS', '18');     % readable when the figure fills a slide
S = plot_style();
Q = load_results(S);
h = Q.h;
% Hysing et al. (2009) TC1, group 1 (TP2D), finest level, from featflow.de:
% columns t, area, circularity, centroid y, rise velocity; shape at t = 3
H  = readmatrix(fullfile(fileparts(mfilename('fullpath')), '..', 'reference_hysing', 'c1g1l7.txt'), 'FileType', 'text');
Hs = readmatrix(fullfile(fileparts(mfilename('fullpath')), '..', 'reference_hysing', 'c1g1l7s.txt'), 'FileType', 'text');

fig = S.figure(1400, 900);
% the bubble takes 2/5 of the width, the three time histories 3/5
tl = tiledlayout(fig, 3, 5, 'TileSpacing', 'compact', 'Padding', 'compact');

ax = nexttile(tl, 1, [3 2]); hold(ax, 'on');
A0 = Q.field(Q.A(:,1)); A1 = Q.field(Q.A(:,end));
S.field(ax, Q.Xq, Q.Yq, Q.field(hypot(Q.UX(:,end), Q.UY(:,end))));
colormap(ax, S.cmap); add_colorbar(ax, S, 'Velocity magnitude |u|');
contour(ax, Q.Xq, Q.Yq, A0, [0.5 0.5], '--', 'LineColor', 'w', 'LineWidth', 1.5);
contour(ax, Q.Xq, Q.Yq, A1, [0.5 0.5], '-', 'LineColor', S.ink, 'LineWidth', 2);
plot(ax, Hs(:,1), Hs(:,2), '.', 'Color', 'm', 'MarkerSize', 4);   % Hysing shape
[Xa, Ya] = meshgrid(linspace(0.05, 0.95, 10), linspace(0.05, 1.95, 20));
FU = scatteredInterpolant(Q.x, Q.y, Q.UX(:,end)); FV = scatteredInterpolant(Q.x, Q.y, Q.UY(:,end));
quiver(ax, Xa, Ya, FU(Xa, Ya), FV(Xa, Ya), 1.2, 'Color', 'w', 'LineWidth', 1);
style_field_axes(ax, S, 1, 2); xticks(ax, 0:0.5:1); yticks(ax, 0:0.5:2);
title(ax, sprintf('Bubble and velocity, t = %.1f', Q.t(end)));   % dashed: t = 0; magenta: Hysing, h = 1/320

ref = struct('t', {[0.9213], [1.9041], [3.0]}, 'v', {0.2417, 0.9013, 1.0813});
rows = {8, 'Rise velocity', 1, 5; 9, 'Circularity', 2, 3; 6, 'Centroid (y)', 3, 4};
for r = 1:3
    ax = nexttile(tl, 5*(r - 1) + 3, [1 3]); hold(ax, 'on');
    c = rows{r, 1};
    plot(ax, h(:,1), h(:,c), '-', 'Color', S.blue, 'LineWidth', 2, 'DisplayName', 'Finch, h = 1/80');
    plot(ax, H(:,1), H(:,rows{r, 4}), '--', 'Color', S.orange, 'LineWidth', 2, ...
         'DisplayName', 'Hysing et al. (2009), h = 1/320');
    S.axes(ax); ylabel(ax, rows{r, 2});
    if r == 3, xlabel(ax, 'Time (t)'); end
    if r == 1, legend(ax, 'Location', 'southeast', 'Box', 'off', 'TextColor', S.ink); end
    xlim(ax, [0 h(end,1)]);
end
S.save(fig, 'Figure2_CaseC_RisingBubble.pdf');
fprintf('written: %s\n', S.file('Figure2_CaseC_RisingBubble.pdf'));
