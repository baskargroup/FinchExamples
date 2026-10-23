% Figure 2: the interface normal on the 2D disk, zoomed on a few cells of
% the interface. Orange: exact normal (towards the centre). Coloured: the
% method. Rows pair each finite difference with its finite element
% counterpart: forward | element-wise (one-sided), central | nodal (centred).
% Element-wise normals are drawn from every element that shares a node, so
% the fan of arrows at a node is the jump between elements.
clc; close all; clear;
plot_common;
fd = S.load('2d_fd.csv');
n1 = S.load('2d_nodal.csv');  e1 = S.load('2d_elem.csv');

th  = 30;                                    % polar angle of the window centre
wc  = [xc + R*cosd(th), yc + R*sind(th)];
win = [wc(1) - 3*h, wc(1) + 3*h, wc(2) - 3*h, wc(2) + 3*h];
L   = 0.8*h;                                 % arrow length

panels = { ...
  'FD forward',                    fd(:,1), fd(:,2), fd(:,4), fd(:,5), C.onesided; ...
  'FEM element-wise',              e1(:,3), e1(:,4), e1(:,5), e1(:,6), C.onesided; ...
  'FD central (Mirjalili et al.)', fd(:,1), fd(:,2), fd(:,8), fd(:,9), C.central; ...
  'FEM nodal',                     n1(:,1), n1(:,2), n1(:,4), n1(:,5), C.central};

fig = S.figure(900, 900);
tl  = tiledlayout(fig, 2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
tt  = linspace(0, 360, 721);
for p = 1:size(panels, 1)
    [name, x, y, gx, gy, col] = panels{p, :};
    ax = nexttile(tl); hold(ax, 'on');
    % element grid and the phi = 0.5 circle
    for g = (ceil(win(1)/h):floor(win(2)/h))*h, xline(ax, g, '-', 'Color', S.grid, 'LineWidth', 0.6); end
    for g = (ceil(win(3)/h):floor(win(4)/h))*h, yline(ax, g, '-', 'Color', S.grid, 'LineWidth', 0.6); end
    plot(ax, xc + R*cosd(tt), yc + R*sind(tt), '-', 'Color', S.ink2, 'LineWidth', 1);
    in = x >= win(1) & x <= win(2) & y >= win(3) & y <= win(4) & ~isnan(gx);
    x = x(in); y = y(in); m = sqrt(gx(in).^2 + gy(in).^2);
    nx = gx(in)./m; ny = gy(in)./m;
    ex = -(x - xc)./hypot(x - xc, y - yc); ey = -(y - yc)./hypot(x - xc, y - yc);
    hE = quiver(ax, x, y, L*ex, L*ey, 0, 'Color', C.exact, 'LineWidth', 1.6, ...
                'MaxHeadSize', 0.4, 'DisplayName', 'Exact');
    hN = quiver(ax, x, y, L*nx, L*ny, 0, 'Color', col, 'LineWidth', 1.3, ...
                'MaxHeadSize', 0.4, 'DisplayName', name);
    plot(ax, x, y, '.', 'Color', S.ink, 'MarkerSize', 8);
    S.axes(ax); grid(ax, 'off');
    axis(ax, 'equal'); xlim(ax, win(1:2)); ylim(ax, win(3:4));
    xlabel(ax, 'Position (x)'); ylabel(ax, 'Position (y)');
    title(ax, name);
    legend(ax, [hE hN], 'Location', 'southoutside', 'Orientation', 'horizontal', ...
           'Box', 'off', 'TextColor', S.ink2);
end
S.save(fig, 'Figure2_Normals2D.pdf');
