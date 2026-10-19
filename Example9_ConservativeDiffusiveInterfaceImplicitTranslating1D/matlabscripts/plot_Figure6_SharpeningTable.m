% Figure 6: every run of the sharpening-term study as a table, one panel per
% mesh. Cells give the max error at T times 1000, shaded on a log scale; the
% smallest error in each column is bold. * = finished but phi left
% [-0.01, 1.01]; div = diverged.
clc; close all; clear;
sweep_common;

rows = {"extrap", 1, 'Extrapolated, 1 solve'};
for f = 1:numel(forms)
    for it = iters
        lab = sprintf('%s, %d it.', upper(extractBefore(forms(f), 2)) + extractAfter(forms(f), 1), it);
        if forms(f) == "lagged" && it == 1; lab = 'Lagged, 1 it. (explicit)'; end
        rows(end+1, :) = {forms(f), it, lab}; %#ok<SAGROW>
    end
end
nr = size(rows, 1);  nc = numel(cfls);
lo = -3;  hi = 0;                       % shading range of log10(error)
cmap = S.ramp(256);

fig = S.figure(1300, 600);
tl  = tiledlayout(fig, 1, numel(meshes), 'TileSpacing', 'compact', 'Padding', 'compact');
for m = 1:numel(meshes)
    E = NaN(nr, nc);  ST = strings(nr, nc);
    for i = 1:nr
        s = sweep_series(res, meshes(m), rows{i,1}, rows{i,2}, cfls);
        E(i,:) = s.err;  ST(i,:) = s.status;
    end
    E(ST == "diverged") = NaN;
    ax = nexttile(tl); hold(ax, 'on');
    L = (min(max(log10(E), lo), hi) - lo)/(hi - lo);
    im = imagesc(ax, L, [0 1]);
    im.AlphaData = ~isnan(E);
    colormap(ax, cmap);
    ax.Color = [0.88 0.88 0.86];        % diverged or missing cells
    [~, best] = min(E, [], 1);
    for i = 1:nr
        for j = 1:nc
            if ST(i,j) == "diverged"
                txt = 'div';
            elseif ST(i,j) == "missing"
                txt = '';
            else
                txt = sprintf('%.3g', 1e3*E(i,j));
                if ST(i,j) == "overshoot"; txt = [txt '*']; end %#ok<AGROW>
            end
            tc = S.ink;  if ~isnan(L(i,j)) && L(i,j) > 0.55; tc = 'w'; end
            fw = 'normal'; if best(j) == i && ~isnan(E(i,j)); fw = 'bold'; end
            text(ax, j, i, txt, 'HorizontalAlignment', 'center', 'Color', tc, ...
                 'FontName', S.font, 'FontSize', S.fs - 1, 'FontWeight', fw);
        end
    end
    % thin separators between the forms
    for yb = [1.5 4.5 7.5 10.5]
        plot(ax, [0.5 nc+0.5], [yb yb], '-', 'Color', 'w', 'LineWidth', 2);
    end
    set(ax, 'YDir', 'reverse', 'FontName', S.font, 'FontSize', S.fs, ...
            'XColor', S.ink2, 'YColor', S.ink2, 'TickLength', [0 0], 'Box', 'off');
    xlim(ax, [0.5 nc+0.5]); ylim(ax, [0.5 nr+0.5]);
    xticks(ax, 1:nc); xticklabels(ax, arrayfun(cfllabel, cfls, 'UniformOutput', false));
    yticks(ax, 1:nr);
    if m == 1; yticklabels(ax, rows(:,3)); else; yticklabels(ax, {}); end
    xlabel(ax, 'CFL number u\Deltat/\Deltax', 'Color', S.ink);
    title(ax, sprintf('%d elements: max error at T (\\times10^{-3})', meshes(m)), ...
          'Color', S.ink, 'FontWeight', 'normal', 'FontSize', S.fs + 2);
    ax.Toolbar.Visible = 'off';
end
S.save(fig, 'Figure6_SharpeningTable.pdf');
