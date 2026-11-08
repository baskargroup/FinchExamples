function draw_frame(tl, S, Xq, Yq, Ak, Ck, Xa, Ya, U, V, cmax, t)
% Two panels: gas fraction | dissolved gas with the bubble outline and the
% expansion flow. Shared by Figure 1 and Video 1.
    ax = nexttile(tl); hold(ax, 'on');
    S.field(ax, Xq, Yq, Ak);
    colormap(ax, S.cmap); clim(ax, [0 1]);
    cb = colorbar(ax); cb.Label.String = 'Gas fraction (\alpha)';
    cb.Label.Color = S.ink; cb.Color = S.ink2; cb.TickDirection = 'in';
    style_field_axes(ax, S);
    title(ax, sprintf('t = %.2f', t));

    ax = nexttile(tl); hold(ax, 'on');
    S.field(ax, Xq, Yq, Ck);
    colormap(ax, S.cmap); clim(ax, [0 cmax]);
    cb = colorbar(ax); cb.Label.String = 'c / c_{sat}';
    cb.Label.Color = S.ink; cb.Color = S.ink2; cb.TickDirection = 'in';
    contour(ax, Xq, Yq, Ak, [0.5 0.5], '-', 'LineColor', S.ink, 'LineWidth', 1.5);
    quiver(ax, Xa, Ya, U, V, 1.5, 'Color', 'w', 'LineWidth', 1);
    style_field_axes(ax, S);
    title(ax, 'Dissolved gas and flow');
end
