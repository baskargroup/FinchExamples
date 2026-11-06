% Video 2 (case C): the rising bubble. Phase field and speed with the flow,
% and rise velocity, circularity and centroid in time against Hysing et al.
% (2009) test case 1, drawn up to the current frame.
clc; close all; clear;
setenv('CDI_CASE', 'C');
setenv('CDI_FS', '14');     % readable when the video fills a slide
S = plot_style();
Q = load_results(S);
h = Q.h;
% Hysing et al. (2009) TC1, TP2D finest level: t, area, circularity, centroid y, rise velocity
H = readmatrix(fullfile(fileparts(mfilename('fullpath')), '..', 'reference_hysing', 'c1g1l7.txt'), 'FileType', 'text');
try
    vid = VideoWriter(S.file('Video2_CaseC_RisingBubble.mp4'), 'MPEG-4');
catch
    vid = VideoWriter(S.file('Video2_CaseC_RisingBubble.avi'), 'Motion JPEG AVI');
end
vid.FrameRate = 12; open(vid);
[Xa, Ya] = meshgrid(linspace(0.05, 0.95, 10), linspace(0.05, 1.95, 20));
umax = max(hypot(Q.UX(:), Q.UY(:)));
% benchmark values: max rise velocity, min circularity, final centroid
ref = {[0.9213 0.2417], [1.9041 0.9013], [3.0 1.0813]};
rows = {8, 'Rise velocity', [0 0.27], 5; 9, 'Circularity', [0.88 1.01], 3; 6, 'Centroid (y)', [0.5 1.15], 4};
fig = S.figure(1600, 740);   % wider than the screen gets cut; -r140 for detail
% CDI_FRAMES limits the number of frames, to check the layout quickly.
nk = numel(Q.t);
if ~isempty(getenv('CDI_FRAMES')), nk = min(nk, str2double(getenv('CDI_FRAMES'))); end
for k = 1:nk
    clf(fig);
    % the two fields take 6/10 of the width (three tiles each), the histories 4/10
    tl = tiledlayout(fig, 3, 10, 'TileSpacing', 'compact', 'Padding', 'compact');
    A1 = Q.field(Q.A(:,k));
    ax = nexttile(tl, 1, [3 3]); hold(ax, 'on');
    S.field(ax, Q.Xq, Q.Yq, A1);
    colormap(ax, S.cmap); clim(ax, [0 1]); add_colorbar(ax, S, 'Phase field (\phi)');
    contour(ax, Q.Xq, Q.Yq, A1, [0.5 0.5], '-', 'LineColor', S.ink, 'LineWidth', 1.5);
    style_field_axes(ax, S, 1, 2); xticks(ax, 0:0.5:1); yticks(ax, 0:0.5:2);
    title(ax, sprintf('Phase field, t = %.2f', Q.t(k)));
    ax = nexttile(tl, 4, [3 3]); hold(ax, 'on');
    S.field(ax, Q.Xq, Q.Yq, Q.field(hypot(Q.UX(:,k), Q.UY(:,k))));
    colormap(ax, S.cmap); clim(ax, [0 umax]); add_colorbar(ax, S, 'Velocity magnitude |u|');
    contour(ax, Q.Xq, Q.Yq, A1, [0.5 0.5], '-', 'LineColor', 'w', 'LineWidth', 1.5);
    FU = scatteredInterpolant(Q.x, Q.y, Q.UX(:,k)); FV = scatteredInterpolant(Q.x, Q.y, Q.UY(:,k));
    quiver(ax, Xa, Ya, FU(Xa, Ya), FV(Xa, Ya), 1.2, 'Color', 'w', 'LineWidth', 1);
    style_field_axes(ax, S, 1, 2); xticks(ax, 0:0.5:1); yticks(ax, 0:0.5:2);
    title(ax, 'Velocity');
    m = h(:,1) <= Q.t(k) + 1e-12;
    for r = 1:3
        ax = nexttile(tl, 10*(r - 1) + 7, [1 4]); hold(ax, 'on');
        plot(ax, h(m,1), h(m,rows{r,1}), '-', 'Color', S.blue, 'LineWidth', 2.5, ...
             'DisplayName', 'Finch, h = 1/80');
        plot(ax, H(:,1), H(:,rows{r,4}), '--', 'Color', S.orange, 'LineWidth', 2, ...
             'DisplayName', 'Hysing et al. (2009), h = 1/320');
        S.axes(ax); ylabel(ax, rows{r,2}); xlim(ax, [0 3]); ylim(ax, rows{r,3});
        if r == 1, legend(ax, 'Location', 'southeast', 'Box', 'off', 'TextColor', S.ink); end
        if r == 3, xlabel(ax, 'Time (t)'); end
    end
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
