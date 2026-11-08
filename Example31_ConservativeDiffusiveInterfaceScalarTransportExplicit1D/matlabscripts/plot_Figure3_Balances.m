clc; close all; clear;
% Figure 3: column edges and the gas volume balance vs time
% edges.csv, per snapshot: t, left edge, right edge (alpha = 0.5 crossings)
% history.csv, per step:   t, int alpha, volume created so far (sum dt int theta15 J),
%                          int c + int alpha/theta15
S = plot_style();
E = S.load('edges.csv');
H = S.load('history.csv');
t = H(:,1);

fig = S.figure(800, 700);
tl  = tiledlayout(fig, 3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

ax = nexttile(tl); hold(ax, 'on');
plot(ax, E(:,1), E(:,2), 'o-', 'Color', S.blue, 'LineWidth', 1.5, 'DisplayName', 'Left edge');
plot(ax, E(:,1), E(:,3), 'o-', 'Color', S.orange, 'LineWidth', 1.5, 'DisplayName', 'Right edge');
S.axes(ax); ylabel(ax, 'Edge (\alpha = 0.5)'); ax.XTickLabel = [];
legend(ax, 'Location', 'east', 'Box', 'off', 'TextColor', S.ink2);

ax = nexttile(tl); hold(ax, 'on');
plot(ax, t, H(:,2) - H(1,2), '-', 'Color', S.blue, 'LineWidth', 2, 'DisplayName', 'Gas volume gained');
plot(ax, t(1:1000:end), H(1:1000:end,3), 'o', 'Color', S.orange, 'LineWidth', 1, ...
     'MarkerSize', 5, 'DisplayName', 'Volume created, \int\theta_{15}J');
S.axes(ax); ylabel(ax, 'Volume'); ax.XTickLabel = [];
legend(ax, 'Location', 'northwest', 'Box', 'off', 'TextColor', S.ink2);

ax = nexttile(tl); hold(ax, 'on');
plot(ax, t, (H(:,2) - H(1,2)) - H(:,3), '-', 'Color', S.blue, 'LineWidth', 2);
S.axes(ax); ylabel(ax, 'Volume balance error'); xlabel(ax, 'Time');
S.save(fig, 'Figure3_Balances.pdf');
