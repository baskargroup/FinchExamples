% Figure 2: interface (phi = 0.5) at every snapshot, colored by time
clc; close all; clear;
plot_common;
cmap = S.ramp(nsnap);

fig = S.figure(620, 520);
ax  = axes(fig); hold(ax, 'on');
for k = 1:nsnap
    F.Values = PHI(:,k);
    contour(ax, Xq, Yq, F(Xq, Yq), [0.5 0.5], 'LineColor', cmap(k,:), 'LineWidth', 1.5);
end
style_field_axes(ax, S);
colormap(ax, cmap); clim(ax, [times(1) times(end)]);
cb = colorbar(ax); cb.Label.String = 'Time'; cb.Label.Color = S.ink;
cb.Color = S.ink2; cb.TickDirection = 'in';
S.save(fig, 'Figure2_TimeSeries.pdf');
