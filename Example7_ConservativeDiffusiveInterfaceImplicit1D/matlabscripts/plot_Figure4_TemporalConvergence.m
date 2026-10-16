% Figure 4: temporal convergence of the semi-implicit scheme.
% Error against a fine-dt reference on the same mesh (spatial error cancels),
% for two Picard iterations per step (as derived) and one. Reads
% output/convergence.csv written by convergence.jl.
clc; close all; clear;
S = plot_style();
d    = S.load('convergence.csv');      % columns: CFL, dt, err(2 iter), err(1 iter)
dt   = d(:,2);  e2 = d(:,3);  e1 = d(:,4);

fig = S.figure(620, 460);
ax  = axes(fig);  hold(ax, 'on');
set(ax, 'XScale', 'log', 'YScale', 'log');

h2 = plot(ax, dt, e2, '-o', 'Color', S.blue,   'LineWidth', 1.4, 'MarkerSize', 6, 'MarkerFaceColor', 'w');
h1 = plot(ax, dt, e1, '-s', 'Color', S.orange, 'LineWidth', 1.4, 'MarkerSize', 6, 'MarkerFaceColor', 'w');

% Slope triangles under each series, spanning one halving of dt in the
% asymptotic range. The horizontal leg is 1 decade-ratio in dt, the vertical
% leg is the expected order.
slope_triangle(ax, dt(4), dt(3), 0.40*e2(4), 2, S);
slope_triangle(ax, dt(4), dt(3), 0.40*e1(4), 1, S);

S.axes(ax);
grid(ax, 'minor');
xlabel(ax, 'Time step \Deltat (\times10^{-3})');
ylabel(ax, 'Max nodal error');
xlim(ax, [dt(end)/1.4, dt(1)*1.4]);
ylim(ax, [min(e2)/4, max(e1)*3]);
xticks(ax, flipud(dt));
xticklabels(ax, arrayfun(@(v) sprintf('%g', 1e3*v), flipud(dt), 'UniformOutput', false));
ax.XMinorTick = 'off';
legend(ax, [h2 h1], {'Two Picard iterations', 'One Picard iteration'}, ...
       'Location', 'northwest', 'Box', 'off', 'FontName', S.font, 'FontSize', S.fs);

% observed orders from least-squares fits, printed for the record
p2 = polyfit(log(dt), log(e2), 1);  p1 = polyfit(log(dt), log(e1), 1);
fprintf('fitted order: two iterations %.2f, one iteration %.2f\n', p2(1), p1(1));

S.save(fig, 'Figure4_TemporalConvergence.pdf');

% Right triangle from (x1,y1) to (x2, y1*(x2/x1)^p) with the legs labelled:
% "1" under the horizontal leg, the order p beside the vertical leg.
function slope_triangle(ax, x1, x2, y1, p, S)
    y2 = y1*(x2/x1)^p;
    c  = S.ink2;
    plot(ax, [x1 x2 x2 x1], [y1 y1 y2 y1], '-', 'Color', c, 'LineWidth', 1.0);
    text(ax, sqrt(x1*x2), y1/1.25, '1', 'Color', c, 'FontName', S.font, 'FontSize', S.fs-1, ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
    text(ax, x2*1.08, sqrt(y1*y2), sprintf('%d', p), 'Color', c, 'FontName', S.font, 'FontSize', S.fs-1, ...
         'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
end
