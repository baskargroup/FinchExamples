% Figure 6 (case C): rising bubble, Finch against Hysing et al. (2009) test
% case 1 at comparable mesh sizes. Their TP2D levels 4-7 are h = 1/40, 1/80,
% 1/160, 1/320 (files c1g1l<n>.txt, shapes c1g1l<n>s.txt). Two figures:
%   a: Finch h = 1/80  against Hysing h = 1/80, with their h = 1/320 for reference
%   b: Finch h = 1/120 against Hysing h = 1/80 and 1/160 (no equal level), and 1/320
% Prints the benchmark quantities of every curve.
clc; close all; clear;
setenv('CDI_CASE', 'C');
setenv('CDI_FS', '18');     % readable when the figure fills a slide
S = plot_style();
ref = fileparts(mfilename('fullpath'));
hy = @(lev) readmatrix(fullfile(ref, '..', 'reference_hysing', sprintf('c1g1l%d.txt', lev)), 'FileType', 'text');
hys = @(lev) readmatrix(fullfile(ref, '..', 'reference_hysing', sprintf('c1g1l%ds.txt', lev)), 'FileType', 'text');
hlev = containers.Map({5, 6, 7}, {'1/80', '1/160', '1/320'});
F80  = load_results(S, 'output_caseC');
F120 = load_results(S, 'output_caseC_120');
grey = [0.45 0.45 0.45];
% the Hysing levels coincide, so each gets its own line style and marker
hsty = containers.Map({5, 6, 7}, {'--', ':', '-.'});
hmrk = containers.Map({5, 6, 7}, {'o', 's', '^'});

figs = {
  'a', {F80},  'Finch, h = 1/80',  [5 7], {S.orange, grey};
  'b', {F120}, 'Finch, h = 1/120', [5 6 7], {S.orange, [160 60 200]/255, grey};
};
for f = 1:size(figs, 1)
    tag = figs{f,1}; Q = figs{f,2}{1}; flab = figs{f,3}; levs = figs{f,4}; hcol = figs{f,5};
    % print extrema
    h = Q.h; [vmax, i] = max(h(:,8)); [cmin, j] = min(h(:,9));
    fprintf('%s: max rise velocity %.4f (t = %.2f), min circularity %.4f (t = %.2f), centroid %.4f\n', flab, vmax, h(i,1), cmin, h(j,1), h(end,6));
    for L = levs
        H = hy(L); [vmax, i] = max(H(:,5)); [cmin, j] = min(H(:,3));
        fprintf('Hysing h = %s: max rise velocity %.4f (t = %.2f), min circularity %.4f (t = %.2f), centroid %.4f\n', hlev(L), vmax, H(i,1), cmin, H(j,1), H(end,4));
    end

    fig = S.figure(1400, 700);
    tl = tiledlayout(fig, 3, 5, 'TileSpacing', 'compact', 'Padding', 'compact');
    ax = nexttile(tl, 1, [3 2]); hold(ax, 'on');
    A1 = Q.field(Q.A(:,end));
    [~, hc] = contour(ax, Q.Xq, Q.Yq, A1, [0.5 0.5], '-', 'LineColor', S.blue, 'LineWidth', 2.5);
    hc.DisplayName = flab;
    for k = 1:numel(levs)
        Hs = hys(levs(k));
        [th, ~] = cart2pol(Hs(:,1) - mean(Hs(:,1)), Hs(:,2) - mean(Hs(:,2))); [~, o] = sort(th);
        L = levs(k); m = numel(o) + 1;
        plot(ax, Hs([o; o(1)],1), Hs([o; o(1)],2), hsty(L), 'Color', hcol{k}, 'LineWidth', 2, ...
             'Marker', hmrk(L), 'MarkerSize', 7, 'MarkerIndices', round(k*m/24):round(m/8):m, ...
             'DisplayName', sprintf('Hysing et al., h = %s', hlev(L)));
    end
    S.axes(ax); xlabel(ax, 'Position (x)'); ylabel(ax, 'Position (y)');
    axis(ax, 'equal'); xlim(ax, [0 1]); ylim(ax, [0.6 1.5]); xticks(ax, 0:0.5:1); yticks(ax, 0.6:0.3:1.5);
    legend(ax, 'Location', 'south', 'Box', 'off', 'TextColor', S.ink);
    title(ax, 'Bubble shape, t = 3.0');

    rows = {8, 'Rise velocity', 5; 9, 'Circularity', 3; 6, 'Centroid (y)', 4};
    for k = 1:3
        ax = nexttile(tl, 5*(k - 1) + 3, [1 3]); hold(ax, 'on');
        h = Q.h;
        plot(ax, h(:,1), h(:,rows{k,1}), '-', 'Color', S.blue, 'LineWidth', 2.5, 'DisplayName', flab);
        for n = 1:numel(levs)
            H = hy(levs(n));
            L = levs(n); m = size(H, 1);
            plot(ax, H(:,1), H(:,rows{k,3}), hsty(L), 'Color', hcol{n}, 'LineWidth', 2, ...
                 'Marker', hmrk(L), 'MarkerSize', 7, 'MarkerIndices', round(n*m/30):round(m/10):m, ...
                 'DisplayName', sprintf('Hysing et al., h = %s', hlev(L)));
        end
        S.axes(ax); ylabel(ax, rows{k,2}); xlim(ax, [0 3]);
        if k == 1, legend(ax, 'Location', 'southeast', 'Box', 'off', 'TextColor', S.ink); end
        if k == 3, xlabel(ax, 'Time (t)'); end
    end
    S.save(fig, sprintf('Figure6%s_CaseC_SameMesh.pdf', tag));
    fprintf('written: %s\n', S.file(sprintf('Figure6%s_CaseC_SameMesh.pdf', tag)));
end
