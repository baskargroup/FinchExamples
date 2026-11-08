function style_field_axes(ax, S)
% Axis styling shared by the 2D field plots.
    S.axes(ax);
    xlabel(ax, 'Position (x)'); ylabel(ax, 'Position (y)');
    axis(ax, 'equal'); xlim(ax, [0 1]); ylim(ax, [0 1]);
    xticks(ax, 0:0.2:1); yticks(ax, 0:0.2:1);
end
