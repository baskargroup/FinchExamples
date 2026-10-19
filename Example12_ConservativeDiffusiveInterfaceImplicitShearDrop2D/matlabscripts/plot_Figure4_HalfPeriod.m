% Figure 4: the most stretched state, t = T/2, with the element grid.
% Left: the whole spiral. Right: zoom on the thin tail, where the filament
% is only a few elements wide.
clc; close all; clear;
plot_common;
data = S.load('solution.csv');
F.Values = data(:,5);  Zh = F(Xq, Yq);

fig = S.figure(1000, 500);
tl  = tiledlayout(fig, 1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

% tail window: around the phi = 0.5 point farthest from the spiral centre
in  = data(:,5) >= 0.5;
xs  = data(in,1); ys = data(in,2);
[~, k] = max((xs - 0.5).^2 + (ys - 0.5).^2);
win = {[0 1 0 1], [xs(k)-0.06 xs(k)+0.06 ys(k)-0.06 ys(k)+0.06]};
ttl = {sprintf('t = T/2 = %.1f', Tp/2), 'The tail of the filament'};
for p = 1:2
    ax = nexttile(tl); hold(ax, 'on');
    contourf(ax, Xq, Yq, Zh, 20, 'LineColor', 'none');
    colormap(ax, S.ramp(64)); clim(ax, [0 1]);
    contour(ax, Xq, Yq, Zh, [0.5 0.5], '-', 'LineColor', S.ink, 'LineWidth', 1.5);
    S.axes(ax);
    axis(ax, 'equal'); xlim(ax, win{p}(1:2)); ylim(ax, win{p}(3:4));
    xlabel(ax, 'Position (x)'); ylabel(ax, 'Position (y)');
    title(ax, ttl{p});
    if p == 1
        rectangle(ax, 'Position', [win{2}(1) win{2}(3) 0.12 0.12], ...
            'EdgeColor', S.orange, 'LineWidth', 1.5);
    else
        % element grid lines
        gx = (floor(win{p}(1)/dx):ceil(win{p}(2)/dx))*dx;
        gy = (floor(win{p}(3)/dx):ceil(win{p}(4)/dx))*dx;
        for g = gx, xline(ax, g, '-', 'Color', [S.grid 0.5], 'LineWidth', 0.4); end
        for g = gy, yline(ax, g, '-', 'Color', [S.grid 0.5], 'LineWidth', 0.4); end
    end
end
cb = colorbar(ax); cb.Layout.Tile = 'east';
cb.Label.String = 'Phase field (\phi)'; cb.Label.Color = S.ink;
cb.Color = S.ink2; cb.TickDirection = 'in';
S.save(fig, 'Figure4_HalfPeriod.pdf');
fprintf('phi at t = T/2: min %.3e, max %.6f\n', min(data(:,5)), max(data(:,5)));
