% Figure 5: where the error is. Left: |phi - exact| at t = T over the domain.
% Right: the phi = 0.5 contour, numerical and exact, zoomed on the leading edge
% with the mesh lines drawn, so sub-cell deviations are visible.
clc; close all; clear;
S = plot_style();
d = S.load('solution.csv');            % x, y, phi, exact
x = d(:,1); y = d(:,2); p = d(:,3); e = d(:,4);
n = round(sqrt(numel(x)));  h = 1/(n-1);
[~, ix] = sort(round(y/h)*n + round(x/h));
X = reshape(x(ix), n, n)'; Y = reshape(y(ix), n, n)'; P = reshape(p(ix), n, n)'; E = reshape(e(ix), n, n)';
xc = 0.7; yc = 0.5; R = 0.15;

fig = S.figure(1000, 440);
tl = tiledlayout(fig, 1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

ax1 = nexttile(tl);
imagesc(ax1, [0 1], [0 1], abs(P - E)); set(ax1, 'YDir', 'normal'); axis(ax1, 'equal', 'tight');
colormap(ax1, flipud(gray)); cb = colorbar(ax1); cb.Label.String = '|\phi - \phi_{exact}|';
hold(ax1, 'on'); th = linspace(0, 2*pi, 400);
plot(ax1, xc + R*cos(th), yc + R*sin(th), '--', 'Color', S.orange, 'LineWidth', 1);
rectangle(ax1, 'Position', [0.83 0.44 0.04 0.12], 'EdgeColor', S.blue, 'LineWidth', 1.2);
S.axes(ax1); grid(ax1, 'off');
xlabel(ax1, 'Position (x)'); ylabel(ax1, 'Position (y)');
title(ax1, sprintf('Error field, max %.3f', max(abs(P(:) - E(:)))));

ax2 = nexttile(tl); hold(ax2, 'on');
xl = [0.83 0.87]; yl = [0.44 0.56];
% mesh lines in the window
for gx = ceil(xl(1)/h)*h : h : xl(2), plot(ax2, [gx gx], yl, '-', 'Color', [0.88 0.88 0.88], 'LineWidth', 0.4); end
for gy = ceil(yl(1)/h)*h : h : yl(2), plot(ax2, xl, [gy gy], '-', 'Color', [0.88 0.88 0.88], 'LineWidth', 0.4); end
[~, hn] = contour(ax2, X, Y, P, [0.5 0.5], 'LineColor', S.blue,   'LineWidth', 1.6);
[~, he] = contour(ax2, X, Y, E, [0.5 0.5], 'LineColor', S.orange, 'LineWidth', 1.6, 'LineStyle', '--');
S.axes(ax2); grid(ax2, 'off'); axis(ax2, 'equal'); xlim(ax2, xl); ylim(ax2, yl);
xlabel(ax2, 'Position (x)'); ylabel(ax2, 'Position (y)');
title(ax2, 'Leading edge, \phi = 0.5 contour (grid lines = elements)');
legend(ax2, [hn he], {'Numerical', 'Exact'}, 'Location', 'southwest', 'Box', 'off', 'FontName', S.font, 'FontSize', S.fs);

S.save(fig, 'Figure5_FrontZoom.pdf');
