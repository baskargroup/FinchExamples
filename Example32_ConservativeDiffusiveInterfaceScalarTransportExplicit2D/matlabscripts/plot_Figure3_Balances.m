% Figure 3: gas volume and gas mole balances vs time
% history.csv, per step: t, int alpha, volume created so far, int c + int alpha/theta15, ...
clc; close all; clear;
S = plot_style();
H = S.load('history.csv');
t = H(:,1);

fig = S.figure(800, 600);
tl  = tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = nexttile(tl); hold(ax, 'on');
plot(ax, t, (H(:,2) - H(1,2)) - H(:,3), '-', 'Color', S.blue, 'LineWidth', 2);
S.axes(ax); ylabel(ax, 'Volume balance error'); ax.XTickLabel = [];
ax = nexttile(tl); hold(ax, 'on');
plot(ax, t, H(:,4) - H(1,4), '-', 'Color', S.orange, 'LineWidth', 2);
S.axes(ax); ylabel(ax, 'Change in gas moles'); xlabel(ax, 'Time');
S.save(fig, 'Figure3_Balances.pdf');
