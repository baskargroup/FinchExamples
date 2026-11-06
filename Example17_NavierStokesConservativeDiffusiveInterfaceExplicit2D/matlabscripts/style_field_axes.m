function style_field_axes(ax, S, Lx, Ly)
% Axis styling shared by the 2D field plots.
    S.axes(ax);
    xlabel(ax, 'Position (x)'); ylabel(ax, 'Position (y)');
    axis(ax, 'equal'); xlim(ax, [0 Lx]); ylim(ax, [0 Ly]);
    xticks(ax, 0:0.25:1); yticks(ax, 0:0.25:2);
    ax.Layer = 'top';
end
