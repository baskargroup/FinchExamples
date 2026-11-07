% Figure 7 (case C): the 80 x 160 mesh of the rising bubble with the initial
% bubble (phi = 0.5 contour), and a zoom on the interface showing the phase
% field across the cells (eps = dx: the tanh profile spans about four cells).
clc; close all; clear;
setenv('CDI_CASE', 'C');
setenv('CDI_FS', '18');
S = plot_style();
Q = load_results(S);
nel = 80; Lx = 1; Ly = 2; dx = Lx / nel;
xc = 0.5; yc = 0.5; R = 0.25;

fig = S.figure(1150, 900);
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% whole mesh
ax = nexttile(tl); hold(ax, 'on');
xs = 0:dx:Lx; ys = 0:dx:Ly;
plot(ax, [xs; xs], [zeros(size(xs)); Ly*ones(size(xs))], '-', 'Color', [0.75 0.75 0.75], 'LineWidth', 0.3);
plot(ax, [zeros(size(ys)); Lx*ones(size(ys))], [ys; ys], '-', 'Color', [0.75 0.75 0.75], 'LineWidth', 0.3);
th = linspace(0, 2*pi, 400);
plot(ax, xc + R*cos(th), yc + R*sin(th), '-', 'Color', S.orange, 'LineWidth', 2.5);
zb = [xc - 0.1, yc + R - 0.1, 0.2, 0.2];            % zoom box on the top of the bubble
rectangle(ax, 'Position', zb, 'EdgeColor', S.blue, 'LineWidth', 2);
S.axes(ax); grid(ax, 'off'); axis(ax, 'equal'); xlim(ax, [0 Lx]); ylim(ax, [0 Ly]);
xticks(ax, 0:0.5:1); yticks(ax, 0:0.5:2);
xlabel(ax, 'Position (x)'); ylabel(ax, 'Position (y)');
title(ax, sprintf('%d x %d elements, \\Deltax = 1/%d', nel, 2*nel, nel));

% zoom: phase field on the cells, interface contour
ax = nexttile(tl); hold(ax, 'on');
A0 = Q.field(Q.A(:,1));
S.field(ax, Q.Xq, Q.Yq, A0); colormap(ax, S.cmap); clim(ax, [0 1]);
add_colorbar(ax, S, 'Phase field (\phi)');
xs = zb(1):dx:zb(1)+zb(3); ys = zb(2):dx:zb(2)+zb(4);
plot(ax, [xs; xs], [zb(2)*ones(size(xs)); (zb(2)+zb(4))*ones(size(xs))], '-', 'Color', 'w', 'LineWidth', 0.8);
plot(ax, [zb(1)*ones(size(ys)); (zb(1)+zb(3))*ones(size(ys))], [ys; ys], '-', 'Color', 'w', 'LineWidth', 0.8);
plot(ax, xc + R*cos(th), yc + R*sin(th), '-', 'Color', 'w', 'LineWidth', 2.5);
S.axes(ax); grid(ax, 'off'); axis(ax, 'equal'); ax.Layer = 'top';
xlim(ax, [zb(1), zb(1)+zb(3)]); ylim(ax, [zb(2), zb(2)+zb(4)]);
xticks(ax, zb(1):0.1:zb(1)+zb(3)); yticks(ax, zb(2):0.1:zb(2)+zb(4));
xlabel(ax, 'Position (x)'); ylabel(ax, 'Position (y)');
title(ax, 'Interface at t = 0, \epsilon = \Deltax');
S.save(fig, 'Figure7_CaseC_Mesh.pdf');
fprintf('written: %s\n', S.file('Figure7_CaseC_Mesh.pdf'));
