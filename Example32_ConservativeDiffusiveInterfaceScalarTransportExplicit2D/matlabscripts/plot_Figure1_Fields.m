% Figure 1: gas fraction and concentration at the final time, with the
% alpha = 0.5 contours of the initial (dashed) and final (solid) bubbles
clc; close all; clear;
plot_common;
A0 = field(A(:,1)); Af = field(A(:,end)); Cf = field(C(:,end));

fig = S.figure(1200, 520);
tl  = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
panels = {Af, Cf}; labels = {'Gas fraction (\alpha)', 'c / c_{sat}'};
for k = 1:2
    ax = nexttile(tl); hold(ax, 'on');
    S.field(ax, Xq, Yq, panels{k});
    colormap(ax, S.cmap);
    if k == 1, clim(ax, [0 1]); end
    cb = colorbar(ax); cb.Label.String = labels{k};
    cb.Label.Color = S.ink; cb.Color = S.ink2; cb.TickDirection = 'in';
    [~, h0] = contour(ax, Xq, Yq, A0, [0.5 0.5], '--', 'LineColor', 'w', ...
                      'LineWidth', 1.5, 'DisplayName', 'Initial');
    [~, hf] = contour(ax, Xq, Yq, Af, [0.5 0.5], '-', 'LineColor', S.ink, ...
                      'LineWidth', 1.5, 'DisplayName', sprintf('t = %g', times(end)));
    style_field_axes(ax, S);
    if k == 1
        legend(ax, [h0 hf], 'Location', 'southeast', 'Box', 'off', 'TextColor', 'w');   % on the dark-blue liquid
    end
end
S.save(fig, 'Figure1_Fields.pdf');
