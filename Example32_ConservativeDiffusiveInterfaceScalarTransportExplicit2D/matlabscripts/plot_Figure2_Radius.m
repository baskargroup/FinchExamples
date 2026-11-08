% Figure 2: equivalent radius sqrt(V/pi) of each bubble vs time
% history.csv, per step: t, int alpha, volume created so far,
% int c + int alpha/theta15, then the gas volume of bubbles A..E
clc; close all; clear;
S = plot_style();
H = S.load('history.csv');
t = H(:,1);
names  = {'A, r_0 = 0.08', 'B, r_0 = 0.06', 'C, r_0 = 0.06, pair', ...
          'D, r_0 = 0.06, pair', 'E, r_0 = 0.04'};
styles = {'-', '-', '--', ':', '-'};
cols   = {S.ink, S.blue, S.orange, S.orange, S.ink2};

fig = S.figure(800, 450);
ax  = axes(fig); hold(ax, 'on');
for k = 1:5
    plot(ax, t, sqrt(H(:,4+k)/pi), styles{k}, 'Color', cols{k}, 'LineWidth', 2, ...
         'DisplayName', names{k});
end
S.axes(ax); xlabel(ax, 'Time'); ylabel(ax, 'Equivalent radius');
legend(ax, 'Location', 'northwest', 'Box', 'off', 'TextColor', S.ink2);
S.save(fig, 'Figure2_Radius.pdf');
