% Figure 3: angle between each method's normal and the exact normal, at every
% node in the interface band (0.05 < phi < 0.95), against the polar angle
% around the disk. Each row pairs a finite difference with its finite element
% counterpart. Element-wise has one point per element sharing the node, so
% its vertical spread at a node is the jump between elements.
clc; close all; clear;
plot_common;
fd = S.load('2d_fd.csv');
n1 = S.load('2d_nodal.csv');  e1 = S.load('2d_elem.csv');

% phi at the element rows, from the nodal file
key = @(x, y) round(x/h) * 1e4 + round(y/h);
pmap = containers.Map(key(n1(:,1), n1(:,2)), n1(:,3));
pe1 = cell2mat(values(pmap, num2cell(key(e1(:,3), e1(:,4)))));

% rows: one-sided FD | element-wise FEM, centred FD | nodal FEM
M = { ...
  'FD forward and backward',       [fd(:,1); fd(:,1)], [fd(:,2); fd(:,2)], [fd(:,3); fd(:,3)], ...
                                   [fd(:,4); fd(:,6)], [fd(:,5); fd(:,7)], C.onesided, 'd'; ...
  'FEM element-wise',              e1(:,3), e1(:,4), pe1,     e1(:,5), e1(:,6), C.onesided, '.'; ...
  'FD central (Mirjalili et al.)', fd(:,1), fd(:,2), fd(:,3), fd(:,8), fd(:,9), C.central, 'o'; ...
  'FEM nodal',                     n1(:,1), n1(:,2), n1(:,3), n1(:,4), n1(:,5), C.central, 'o'};

fig = S.figure(1150, 600);
tl  = tiledlayout(fig, 2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
fprintf('%-32s %6s %10s %10s\n', 'method', 'n', 'max (deg)', 'rms (deg)');
for k = 1:size(M, 1)
    [name, x, y, p, gx, gy, col, mk] = M{k, :};
    ok = inband(p) & ~isnan(gx);
    a  = angerr(gx(ok), gy(ok), x(ok), y(ok));
    th = polar(x(ok), y(ok));
    fprintf('%-32s %6d %10.2f %10.2f\n', name, numel(a), max(abs(a)), sqrt(mean(a.^2)));
    ax = nexttile(tl); hold(ax, 'on');
    yline(ax, 0, '-', 'Color', C.exact, 'LineWidth', 1.5);
    plot(ax, th, a, mk, 'Color', col, 'MarkerSize', 4 + 8*strcmp(mk, '.'), 'LineWidth', 0.8);
    S.axes(ax);
    xlim(ax, [0 360]); ylim(ax, [-22 22]);
    xticks(ax, 0:45:360); yticks(ax, -20:10:20);
    xlabel(ax, 'Polar angle (deg)'); ylabel(ax, 'Angle error (deg)');
    title(ax, sprintf('%s:  max %.2f^\\circ, rms %.2f^\\circ', name, max(abs(a)), sqrt(mean(a.^2))));
end
S.save(fig, 'Figure3_AngleError.pdf');
