function style_field_axes(ax, S)
% Axis styling shared by the 2D field plots (Figure 1, Figure 2, Video 1).
    S.axes(ax);
    xlabel(ax, 'Position (x)'); ylabel(ax, 'Position (y)');
    axis(ax, 'equal'); xlim(ax, [0 1]); ylim(ax, [0 1]);
    xticks(ax, 0:0.2:1); yticks(ax, 0:0.2:1);
    grid(ax, 'minor');
    ax.MinorGridColor = S.grid; ax.MinorGridAlpha = 0.6; ax.MinorGridLineStyle = ':';
    ax.XAxis.MinorTickValues = 0:0.05:1;
    ax.YAxis.MinorTickValues = 0:0.05:1;
end
