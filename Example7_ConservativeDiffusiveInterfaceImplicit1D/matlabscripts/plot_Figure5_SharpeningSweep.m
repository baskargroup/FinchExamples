% Figure 5: max error at T against CFL for every treatment of the sharpening
% term, one panel per mesh and number of Picard iterations. The extrapolated
% form takes one solve per step and is repeated in every panel for reference.
% Open markers: the run finished but phi left [-0.01, 1.01]. Crosses on the
% top edge: the run diverged.
clc; close all; clear;
sweep_common;
ytop = 1.5;

fig = S.figure(1250, 760);
tl  = tiledlayout(fig, numel(meshes), numel(iters), 'TileSpacing', 'compact', 'Padding', 'compact');
for m = 1:numel(meshes)
    for it = iters
        ax = nexttile(tl); hold(ax, 'on');
        set(ax, 'XScale', 'log', 'YScale', 'log');
        h = gobjects(0);
        h(end+1) = plot_series(ax, sweep_series(res, meshes(m), "extrap", 1, cfls), ...
                               cfls, xcolor, '--', 'none', ytop);
        for f = 1:numel(forms)
            h(end+1) = plot_series(ax, sweep_series(res, meshes(m), forms(f), it, cfls), ...
                                   cfls, fcolors(f,:), '-', fmarks(f), ytop); %#ok<SAGROW>
        end
        S.axes(ax);
        xlim(ax, [min(cfls)/1.3, max(cfls)*1.3]); ylim(ax, [1e-3 ytop]);
        xticks(ax, cfls); xticklabels(ax, arrayfun(cfllabel, cfls, 'UniformOutput', false));
        ax.XMinorTick = 'off'; ax.XMinorGrid = 'off';
        if it == 1
            itxt = '1 iteration (lagged = explicit \phi^o(1-\phi^o))';
        else
            itxt = sprintf('%d Picard iterations', it);
        end
        title(ax, sprintf('%d elements, %s', meshes(m), itxt));
        if m == numel(meshes); xlabel(ax, 'CFL number u\Deltat/\Deltax'); end
        if it == 1; ylabel(ax, 'Max nodal error at T'); end
        if m == 1 && it == 1
            lg = legend(ax, h, [xname, fnames], 'Box', 'off', 'TextColor', S.ink2, ...
                        'FontName', S.font, 'FontSize', S.fs - 1, 'NumColumns', 5);
            lg.Layout.Tile = 'north';
        end
    end
end
S.save(fig, 'Figure5_SharpeningSweep.pdf');

% One series: filled markers for bounded runs, open for overshoot, a cross on
% the top edge for diverged runs. Returns the line handle for the legend.
function h = plot_series(ax, s, cfls, c, ls, mk, ytop)
    e = s.err; e(s.status == "diverged") = NaN;
    h = plot(ax, cfls, e, ls, 'Color', c, 'LineWidth', 1.6, 'Marker', mk, ...
             'MarkerSize', 7, 'MarkerFaceColor', c, 'MarkerEdgeColor', c);
    over = s.status == "overshoot";
    if any(over) && mk ~= "none"
        plot(ax, cfls(over), e(over), mk, 'MarkerSize', 7, 'MarkerFaceColor', 'w', ...
             'MarkerEdgeColor', c, 'LineWidth', 1.4, 'HandleVisibility', 'off');
    end
    div = s.status == "diverged";
    if any(div)
        plot(ax, cfls(div), ytop/1.15*ones(1, nnz(div)), 'x', 'Color', c, ...
             'MarkerSize', 10, 'LineWidth', 2, 'HandleVisibility', 'off');
    end
end
