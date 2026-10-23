% Figure 1: d(phi)/dx of the 1D tanh profile near the interface.
% Left: finite differences at the nodes. Right: linear finite elements,
% element-wise (one constant segment per element) and nodal. Element-wise =
% forward/backward difference of that element, so the two ends of the
% segments meeting at a node are the forward and backward markers on the
% left; nodal = central difference (main.jl checks both to round-off).
clc; close all; clear;
plot_common;
fd = S.load('1d_fd.csv');
n1 = S.load('1d_nodal.csv');  e1 = S.load('1d_elem.csv');
xs = linspace(x0 - 5*h, x0 + 5*h, 600);
ex = -0.25/eps*sech((xs - x0)/(2*eps)).^2;
win = [x0 - 5*h, x0 + 5*h];

fig = S.figure(1100, 480);
tl  = tiledlayout(fig, 1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

% finite differences
ax = nexttile(tl); hold(ax, 'on');
hE = plot(ax, xs, ex, '--', 'Color', C.exact, 'LineWidth', 2, 'DisplayName', 'Exact');
hF = plot(ax, fd(:,1), fd(:,3), '>', 'Color', C.onesided, 'MarkerSize', 8, ...
          'LineWidth', 1.2, 'DisplayName', 'Forward  (\phi_{i+1} - \phi_i)/h');
hB = plot(ax, fd(:,1), fd(:,4), '<', 'Color', C.onesided, 'MarkerSize', 8, ...
          'LineWidth', 1.2, 'DisplayName', 'Backward  (\phi_i - \phi_{i-1})/h');
hC = plot(ax, fd(:,1), fd(:,5), 'o', 'Color', C.central, 'MarkerSize', 8, ...
          'LineWidth', 1.5, 'DisplayName', 'Central  (\phi_{i+1} - \phi_{i-1})/2h');
decorate(ax, S, win, h, 'Finite differences (at the nodes)');
legend(ax, [hE hF hB hC], 'Location', 'southwest', 'Box', 'off', 'TextColor', S.ink2);

% finite elements
ax = nexttile(tl); hold(ax, 'on');
hE = plot(ax, xs, ex, '--', 'Color', C.exact, 'LineWidth', 2, 'DisplayName', 'Exact');
h1 = segments(ax, e1, C.onesided, 'Element-wise (one value per element)');
hN1 = plot(ax, n1(:,1), n1(:,3), 'o', 'Color', C.central, 'MarkerSize', 8, ...
           'LineWidth', 1.5, 'DisplayName', 'Nodal (L2 projection)');
decorate(ax, S, win, h, 'Finite elements, linear');
legend(ax, [hE h1 hN1], 'Location', 'southwest', 'Box', 'off', 'TextColor', S.ink2);
S.save(fig, 'Figure1_Gradient1D.pdf');

% equivalences, as checked in main.jl
ok = ~isnan(fd(:,5));
fprintf('max |nodal - central| = %.2e\n', max(abs(n1(ok,3) - fd(ok,5))));

function hs = segments(ax, e, col, name)
% one line per element through its nodes' element-wise values
    els = unique(e(:,1));
    for k = 1:numel(els)
        r = e(:,1) == els(k);
        [xx, o] = sort(e(r,3)); g = e(r,4); g = g(o);
        hs = plot(ax, xx, g, '-', 'Color', col, 'LineWidth', 2.5, 'DisplayName', name);
        plot(ax, xx([1 end]), g([1 end]), '.', 'Color', col, 'MarkerSize', 14, ...
             'HandleVisibility', 'off');
    end
end

function decorate(ax, S, win, h, ttl)
    S.axes(ax);
    xlim(ax, win); ylim(ax, [-11 1]);      % exact minimum is -1/(4 eps) = -10
    gx = (ceil(win(1)/h):floor(win(2)/h))*h;
    for g = gx, xline(ax, g, '-', 'Color', [S.grid 0.8], 'LineWidth', 0.5, 'HandleVisibility', 'off'); end
    xlabel(ax, 'Position (x)'); ylabel(ax, 'd\phi/dx');
    title(ax, ttl);
end
