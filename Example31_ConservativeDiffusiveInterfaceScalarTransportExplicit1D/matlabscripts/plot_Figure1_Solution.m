clc; close all; clear;
% Figure 1: final state, gas fraction, concentration and velocity
% solution.csv: x, alpha, c, u, J (in order of x)
S    = plot_style();
data = S.load('solution.csv');
x = data(:,1); a = data(:,2); c = data(:,3); u = data(:,4);
T = 1.0;            % final time (must match main.jl)

fig = S.figure(800, 700);
tl  = tiledlayout(fig, 3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
labels = {'Gas fraction (\alpha)', 'c / c_{sat}', 'Velocity (u)'};
vals   = {a, c, u};
for k = 1:3
    ax = nexttile(tl); hold(ax, 'on');
    plot(ax, x, vals{k}, 'o-', 'Color', S.blue, 'LineWidth', 1, 'MarkerSize', 3);
    S.axes(ax); xlim(ax, [0 1]); xticks(ax, 0:0.1:1);
    ylabel(ax, labels{k});
    if k < 3, ax.XTickLabel = []; end
end
ylim(nexttile(tl, 1), [-0.05 1.05]); ylim(nexttile(tl, 2), [-0.1 3.1]);
xlabel(ax, 'Position (x)');
title(tl, sprintf('t = %g', T), 'FontName', S.font, 'FontSize', S.fs + 2, 'Color', S.ink);
S.save(fig, 'Figure1_Solution.pdf');
