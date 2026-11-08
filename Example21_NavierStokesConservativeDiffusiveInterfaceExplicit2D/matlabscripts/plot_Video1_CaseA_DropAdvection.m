% Video 1 (case A): the drop and the velocity error, consistent (left) and
% inconsistent (right) momentum transport, side by side in time.
clc; close all; clear;
setenv('CDI_CASE', 'A');
setenv('CDI_FS', '14');     % readable when the video fills a slide
S = plot_style();
runs   = {'output_caseA', 'output_caseA_theta23_0'};
labels = {'With consistent term', 'Without consistent term'};
R = cellfun(@(f) load_results(S, f), runs, 'UniformOutput', false);
try
    vid = VideoWriter(S.file('Video1_CaseA_DropAdvection.mp4'), 'MPEG-4');
catch
    vid = VideoWriter(S.file('Video1_CaseA_DropAdvection.avi'), 'Motion JPEG AVI');
end
vid.FrameRate = 10; open(vid);
fig = S.figure(1400, 720);   % wider than the screen gets cut; -r140 for detail
cols = {S.blue, S.orange};
nk = max(numel(R{1}.t), numel(R{2}.t));
% CDI_FRAMES limits the number of frames, to check the layout quickly.
if ~isempty(getenv('CDI_FRAMES')), nk = min(nk, str2double(getenv('CDI_FRAMES'))); end
for k = 1:nk
    clf(fig);
    % fields take 2/3 of the width (two tiles each), the error history 1/3
    tl = tiledlayout(fig, 2, 6, 'TileSpacing', 'compact', 'Padding', 'compact');
    for r = 1:2
        Q = R{r}; kk = min(k, numel(Q.t));
        A1 = Q.field(Q.A(:,kk));
        ax = nexttile(tl, 6*(r - 1) + 1, [1 2]); hold(ax, 'on');
        S.field(ax, Q.Xq, Q.Yq, A1);
        colormap(ax, S.cmap); clim(ax, [0 1]); add_colorbar(ax, S, 'Phase field (\phi)');
        contour(ax, Q.Xq, Q.Yq, A1, [0.5 0.5], '-', 'LineColor', S.ink, 'LineWidth', 1.5);
        style_field_axes(ax, S, Q.Lx, Q.Ly); xticks(ax, 0:0.5:1); yticks(ax, 0:0.5:1);
        title(ax, {labels{r}, sprintf('t = %.2f', Q.t(kk))});
        ax = nexttile(tl, 6*(r - 1) + 3, [1 2]); hold(ax, 'on');
        E = Q.field(log10(max(hypot(Q.UX(:,kk) - 1, Q.UY(:,kk)), 1e-12)));
        S.field(ax, Q.Xq, Q.Yq, E);
        colormap(ax, S.cmap); clim(ax, [-10 1]);
        add_colorbar(ax, S, 'log_{10} |u - (1, 0)|');
        contour(ax, Q.Xq, Q.Yq, A1, [0.5 0.5], '-', 'LineColor', 'w', 'LineWidth', 1.5);
        style_field_axes(ax, S, Q.Lx, Q.Ly); xticks(ax, 0:0.5:1); yticks(ax, 0:0.5:1);
        title(ax, 'Velocity error');
    end
    % error in time, drawn up to this frame
    ax = nexttile(tl, 5, [2 2]); hold(ax, 'on');
    tk = R{1}.t(min(k, numel(R{1}.t)));
    for r = 1:2
        h = R{r}.h; m = h(:,1) <= tk + 1e-12;
        semilogy(ax, h(m,1), max(h(m,10), 1e-16), '-', 'Color', cols{r}, 'LineWidth', 2.5, ...
                 'DisplayName', labels{r});
    end
    set(ax, 'YScale', 'log'); xlim(ax, [0 1]); ylim(ax, [1e-12 1e1]);
    S.axes(ax); xlabel(ax, 'Time (t)'); ylabel(ax, '||u - (1, 0)||_2');
    legend(ax, 'Location', 'southeast', 'Box', 'off', 'TextColor', S.ink);
    title(ax, 'Velocity error (exact: 0)');
    drawnow;
    frame = print(fig, '-RGBImage', '-r140');
    if k == 1, vsz = ceil([size(frame,1) size(frame,2)]/16)*16; end
    frame = frame(1:min(end, vsz(1)), 1:min(end, vsz(2)), :);   % later frames can grow a few pixels
    padded = 255*ones(vsz(1), vsz(2), 3, 'uint8');
    padded(1:size(frame,1), 1:size(frame,2), :) = frame;
    writeVideo(vid, padded);
end
close(vid);
fprintf('video written: %s\n', fullfile(vid.Path, vid.Filename));
