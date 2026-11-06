% Figure 4 (case C): rising bubble on two meshes, 80 x 160 and 120 x 240
% (output_caseC and output_caseC_120), against Hysing et al. (2009) test
% case 1. Left: bubble shapes at t = 3. Right: rise velocity, circularity
% and centroid in time. Prints the benchmark quantities of both runs.
clc; close all; clear;
setenv('CDI_CASE', 'C');
setenv('CDI_FS', '18');     % readable when the figure fills a slide
S = plot_style();
runs   = {'output_caseC', 'output_caseC_120'};
labels = {'Finch, h = 1/80', 'Finch, h = 1/120'};
cols   = {S.blue, [ 60 160  60]/255};
R = cellfun(@(f) load_results(S, f), runs, 'UniformOutput', false);
H  = readmatrix(fullfile(fileparts(mfilename('fullpath')), '..', 'reference_hysing', 'c1g1l7.txt'), 'FileType', 'text');
Hs = readmatrix(fullfile(fileparts(mfilename('fullpath')), '..', 'reference_hysing', 'c1g1l7s.txt'), 'FileType', 'text');

for r = 1:2
    h = R{r}.h;
    [vmax, i] = max(h(:,8)); [cmin, j] = min(h(:,9));
    fprintf('%s: max rise velocity %.4f at t = %.2f, min circularity %.4f at t = %.2f, centroid %.4f at t = %.1f\n', ...
            labels{r}, vmax, h(i,1), cmin, h(j,1), h(end,6), h(end,1));
end

fig = S.figure(1400, 700);
tl = tiledlayout(fig, 3, 5, 'TileSpacing', 'compact', 'Padding', 'compact');

ax = nexttile(tl, 1, [3 2]); hold(ax, 'on');
for r = 1:2
    Q = R{r}; A1 = Q.field(Q.A(:,end));
    [~, hc] = contour(ax, Q.Xq, Q.Yq, A1, [0.5 0.5], '-', 'LineColor', cols{r}, 'LineWidth', 2);
    hc.DisplayName = labels{r};
end
% the benchmark shape is an unordered point set: draw it as a dashed line by angle
[th, ~] = cart2pol(Hs(:,1) - mean(Hs(:,1)), Hs(:,2) - mean(Hs(:,2))); [~, o] = sort(th);
plot(ax, Hs([o; o(1)],1), Hs([o; o(1)],2), '--', 'Color', S.orange, 'LineWidth', 2, ...
     'DisplayName', 'Hysing et al. (2009), h = 1/320');
S.axes(ax); xlabel(ax, 'Position (x)'); ylabel(ax, 'Position (y)');
axis(ax, 'equal'); xlim(ax, [0 1]); ylim(ax, [0.6 1.5]); xticks(ax, 0:0.5:1); yticks(ax, 0.6:0.3:1.5);
legend(ax, 'Location', 'south', 'Box', 'off', 'TextColor', S.ink);
title(ax, 'Bubble shape, t = 3.0');

rows = {8, 'Rise velocity', 5; 9, 'Circularity', 3; 6, 'Centroid (y)', 4};
for k = 1:3
    ax = nexttile(tl, 5*(k - 1) + 3, [1 3]); hold(ax, 'on');
    for r = 1:2
        h = R{r}.h;
        plot(ax, h(:,1), h(:,rows{k,1}), '-', 'Color', cols{r}, 'LineWidth', 2, 'DisplayName', labels{r});
    end
    plot(ax, H(:,1), H(:,rows{k,3}), '--', 'Color', S.orange, 'LineWidth', 2, ...
         'DisplayName', 'Hysing et al. (2009), h = 1/320');
    S.axes(ax); ylabel(ax, rows{k,2}); xlim(ax, [0 3]);
    if k == 1, legend(ax, 'Location', 'southeast', 'Box', 'off', 'TextColor', S.ink); end
    if k == 3, xlabel(ax, 'Time (t)'); end
end
S.save(fig, 'Figure4_CaseC_MeshRefinement.pdf');
fprintf('written: %s\n', S.file('Figure4_CaseC_MeshRefinement.pdf'));
