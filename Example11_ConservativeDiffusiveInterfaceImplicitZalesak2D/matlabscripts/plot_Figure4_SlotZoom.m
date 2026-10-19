% Figure 4: zoom on the slotted disk after one revolution, with the element grid.
% Left: the whole disk. Right: the slot, where the shape changes most.
clc; close all; clear;
plot_common;
F.Values = PHI(:,end);  Zf = F(Xq, Yq);
F.Values = PHI(:,1);    Z0 = F(Xq, Yq);

fig = S.figure(1000, 500);
tl  = tiledlayout(fig, 1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

win = {[xc-0.2 xc+0.2 yc-0.2 yc+0.2], [xc-0.06 xc+0.06 sy1-0.03 sy2+0.03]};
ttl = {'Slotted disk after one revolution', 'The slot'};
for p = 1:2
    ax = nexttile(tl); hold(ax, 'on');
    [~, hN] = contour(ax, Xq, Yq, Zf, [0.5 0.5], '-', 'LineColor', S.blue, ...
        'LineWidth', 1.8, 'DisplayName', 'Numerical');
    [~, hE] = contour(ax, Xq, Yq, Z0, [0.5 0.5], '--', 'LineColor', S.orange, ...
        'LineWidth', 1.8, 'DisplayName', 'Exact (initial field)');
    S.axes(ax);
    axis(ax, 'equal'); xlim(ax, win{p}(1:2)); ylim(ax, win{p}(3:4));
    xlabel(ax, 'Position (x)'); ylabel(ax, 'Position (y)');
    title(ax, ttl{p}, 'FontWeight', 'normal', 'Color', S.ink);
    if p == 2
        % element grid lines
        gx = (floor(win{p}(1)/dx):ceil(win{p}(2)/dx))*dx;
        gy = (floor(win{p}(3)/dx):ceil(win{p}(4)/dx))*dx;
        for g = gx, xline(ax, g, '-', 'Color', [S.grid 0.9], 'LineWidth', 0.4); end
        for g = gy, yline(ax, g, '-', 'Color', [S.grid 0.9], 'LineWidth', 0.4); end
        uistack(hN, 'top'); uistack(hE, 'top');
    end
    legend(ax, [hE hN], 'Location', 'southoutside', 'Orientation', 'horizontal', ...
        'Box', 'off', 'TextColor', S.ink2);
end
S.save(fig, 'Figure4_SlotZoom.pdf');

% Slot width at mid-height, from the phi = 0.5 crossings along y = (sy1+sy2)/2
ym = 0.5*(sy1 + sy2);
xs = linspace(xc - 0.1, xc + 0.1, 2001);
for which = {'initial', 'final'}
    if strcmp(which{1}, 'initial'), F.Values = PHI(:,1); else, F.Values = PHI(:,end); end
    p = F(xs, ym*ones(size(xs)));
    i = find(p(1:end-1) >= 0.5 & p(2:end) < 0.5, 1);
    j = find(p(1:end-1) < 0.5 & p(2:end) >= 0.5, 1, 'last');
    xl = interp1(p([i i+1]), xs([i i+1]), 0.5);
    xr = interp1(p([j j+1]), xs([j j+1]), 0.5);
    fprintf('slot width at mid-height, %s: %.4f (exact %.4f)\n', which{1}, xr - xl, sw);
end
