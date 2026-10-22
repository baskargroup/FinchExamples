% Figure 2: the phase field at every eighth of the period, t = 0, T/8, ..., T.
% The drop is stretched into a spiral up to t = T/2, then the flow reverses.
clc; close all; clear;
plot_common;
keep = round(linspace(1, nsnap, 9));

fig = S.figure(1000, 1000);
tl  = tiledlayout(fig, 3, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
for j = 1:numel(keep)
    ax = nexttile(tl); hold(ax, 'on');
    F.Values = PHI(:,keep(j));
    Z = F(Xq, Yq);
    contourf(ax, Xq, Yq, Z, 20, 'LineColor', 'none');
    colormap(ax, S.ramp(64)); clim(ax, [0 1]);
    contour(ax, Xq, Yq, Z, [0.5 0.5], '-', 'LineColor', S.ink, 'LineWidth', 1.2);
    style_field_axes(ax, S);
    title(ax, sprintf('t = %.2f', times(keep(j))));
end
cb = colorbar(ax); cb.Layout.Tile = 'east';
cb.Label.String = 'Phase field (\phi)'; cb.Label.Color = S.ink;
cb.Color = S.ink2; cb.TickDirection = 'in';
if nx > 400
    % nine filled contours on a fine grid are too much for a vector PDF
    % (the export comes out blank), so write an image PDF instead
    exportgraphics(fig, S.file('Figure2_TimeSeries.pdf'), ...
                   'ContentType', 'image', 'Resolution', 300, 'BackgroundColor', 'w');
else
    S.save(fig, 'Figure2_TimeSeries.pdf');
end
