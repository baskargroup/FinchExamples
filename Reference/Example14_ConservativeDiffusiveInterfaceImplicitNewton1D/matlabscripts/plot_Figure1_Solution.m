clc; close all; clear;
% Figure 1: final state, numerical vs exact translated tanh
% Final state of the 1D CDI solution: numerical vs exact translated tanh.
% solution.csv: x, phi (numerical), phi (exact)
S = plot_style();
data   = S.load('solution.csv');
x      = data(:,1);
phi    = data(:,2);
phi_ex = data(:,3);
T      = 0.25;      % final time (must match main.jl)

fig = S.figure(800, 450);
ax  = axes(fig); hold(ax, 'on');
% markers first, exact line on top so it stays visible through the markers
hN = plot(ax, x, phi, 'o', 'Color', S.blue, 'LineWidth', 1, ...
     'MarkerSize', 5, 'MarkerFaceColor', 'none', ...
     'DisplayName', sprintf('Numerical (t = %g)', T));
hE = plot(ax, x, phi_ex, '-', 'Color', S.orange, 'LineWidth', 2, ...
     'DisplayName', sprintf('Exact (t = %g)', T));
S.axes(ax);
xlabel(ax, 'Position (x)');
ylabel(ax, 'Phase field (\phi)');
xlim(ax, [0 1]); ylim(ax, [-0.05 1.05]);
xticks(ax, 0:0.2:1); yticks(ax, 0:0.25:1);
grid(ax, 'minor');                       % major + minor grid
ax.MinorGridColor = S.grid; ax.MinorGridAlpha = 0.6; ax.MinorGridLineStyle = ':';
ax.XAxis.MinorTickValues = 0:0.05:1;
ax.YAxis.MinorTickValues = 0:0.05:1;
legend(ax, [hE hN], 'Location', 'northwest', 'Box', 'off', 'TextColor', S.ink2);
S.save(fig, 'Figure1_Solution.pdf');
fprintf('max error = %.4e\n', max(abs(phi - phi_ex)));
