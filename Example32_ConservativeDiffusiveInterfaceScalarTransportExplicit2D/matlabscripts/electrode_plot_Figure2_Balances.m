% Figure 2: gas volume against the volume created, and the gas moles against
% what the electrode put in
% history.csv, per step: t, int alpha, volume created so far,
% int c + int alpha/theta15, gas put in by the electrode so far
clc; close all; clear;
S = plot_style('electrode');
H = S.load('history.csv');
t = H(:,1); every = max(1, floor(numel(t)/25));

fig = S.figure(800, 800);
tl  = tiledlayout(fig, 3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = nexttile(tl); hold(ax, 'on');
plot(ax, t, H(:,2) - H(1,2), '-', 'Color', S.blue, 'LineWidth', 2, 'DisplayName', 'Gas volume gained');
plot(ax, t(1:every:end), H(1:every:end,3), 'o', 'Color', S.orange, 'LineWidth', 1, ...
     'MarkerSize', 5, 'DisplayName', 'Volume created, \int\theta_{15}J');
S.axes(ax); ylabel(ax, 'Volume'); ax.XTickLabel = [];
legend(ax, 'Location', 'northwest', 'Box', 'off', 'TextColor', S.ink2);

ax = nexttile(tl); hold(ax, 'on');
plot(ax, t, (H(:,2) - H(1,2)) - H(:,3), '-', 'Color', S.blue, 'LineWidth', 2);
S.axes(ax); ylabel(ax, 'Volume balance error'); ax.XTickLabel = [];

ax = nexttile(tl); hold(ax, 'on');
plot(ax, t, H(:,5), '-', 'Color', S.ink2, 'LineWidth', 2, 'DisplayName', 'Put in by the electrode');
plot(ax, t, H(:,4) - H(1,4), '-', 'Color', S.orange, 'LineWidth', 2, 'DisplayName', 'Gas moles gained');
S.axes(ax); ylabel(ax, 'Gas moles'); xlabel(ax, 'Time');
legend(ax, 'Location', 'northwest', 'Box', 'off', 'TextColor', S.ink2);
S.save(fig, 'Figure2_Balances.pdf');
