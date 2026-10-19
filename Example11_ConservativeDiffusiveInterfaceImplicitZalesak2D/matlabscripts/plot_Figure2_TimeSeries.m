% Figure 2: interface (phi = 0.5) at every eighth of the revolution, colored by time
clc; close all; clear;
plot_common;
keep = round(linspace(1, nsnap, 9));      % t = 0, T/8, ..., T
cmap = S.ramp(numel(keep));

fig = S.figure(620, 520);
ax  = axes(fig); hold(ax, 'on');
for j = 1:numel(keep)
    F.Values = PHI(:,keep(j));
    contour(ax, Xq, Yq, F(Xq, Yq), [0.5 0.5], 'LineColor', cmap(j,:), 'LineWidth', 1.5);
end
style_field_axes(ax, S);
colormap(ax, cmap); clim(ax, [times(1) times(end)]);
cb = colorbar(ax); cb.Label.String = 'Time'; cb.Label.Color = S.ink;
cb.Color = S.ink2; cb.TickDirection = 'in';
S.save(fig, 'Figure2_TimeSeries.pdf');
