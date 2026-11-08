% Video 3 (case D): dam break. Phase field and velocity in time, with the
% front position and column height drawn up to the current frame against
% the Martin and Moyce (1952) front.
clc; close all; clear;
setenv('CDI_CASE', 'D');
setenv('CDI_FS', '14');
S = plot_style();
Q = load_results(S);
a = 0.25; g = 1; Tsc = sqrt(2 * g / a);
F = readmatrix(S.file('front.csv'));          % t, x_front, h_left, T, Z, H
E = readmatrix(fullfile(fileparts(mfilename('fullpath')), '..', 'reference_martin_moyce', 'front_n2.txt'), 'FileType', 'text', 'CommentStyle', '#');
try
    vid = VideoWriter(S.file('Video3_CaseD_DamBreak.mp4'), 'MPEG-4');
catch
    vid = VideoWriter(S.file('Video3_CaseD_DamBreak.avi'), 'Motion JPEG AVI');
end
vid.FrameRate = 10; open(vid);
[Xa, Ya] = meshgrid(linspace(0.03, 0.97, 16), linspace(0.03, 0.72, 12));
umax = max(hypot(Q.UX(:), Q.UY(:)));
fig = S.figure(1400, 760);
% CDI_FRAMES limits the number of frames, to check the layout quickly.
nk = numel(Q.t);
if ~isempty(getenv('CDI_FRAMES')), nk = min(nk, str2double(getenv('CDI_FRAMES'))); end
for k = 1:nk
    clf(fig);
    tl = tiledlayout(fig, 2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
    A1 = Q.field(Q.A(:,k));
    ax = nexttile(tl, 1, [1 2]); hold(ax, 'on');
    S.field(ax, Q.Xq, Q.Yq, A1); colormap(ax, S.cmap); clim(ax, [0 1]); add_colorbar(ax, S, 'Phase field (\phi)');
    contour(ax, Q.Xq, Q.Yq, A1, [0.5 0.5], '-', 'LineColor', 'w', 'LineWidth', 1.5);
    S.axes(ax); axis(ax, 'equal'); xlim(ax, [0 Q.Lx]); ylim(ax, [0 Q.Ly]); ax.Layer = 'top'; grid(ax, 'off');
    xticks(ax, 0:0.25:1); yticks(ax, 0:0.25:0.75); xlabel(ax, 'Position (x)'); ylabel(ax, 'Position (y)');
    title(ax, sprintf('Phase field, T = %.2f', Q.t(k) * Tsc));
    ax = nexttile(tl, 3, [1 2]); hold(ax, 'on');
    S.field(ax, Q.Xq, Q.Yq, Q.field(hypot(Q.UX(:,k), Q.UY(:,k))));
    colormap(ax, S.cmap); clim(ax, [0 umax]); add_colorbar(ax, S, 'Velocity magnitude |u|');
    contour(ax, Q.Xq, Q.Yq, A1, [0.5 0.5], '-', 'LineColor', 'w', 'LineWidth', 1.5);
    FU = scatteredInterpolant(Q.x, Q.y, Q.UX(:,k)); FV = scatteredInterpolant(Q.x, Q.y, Q.UY(:,k));
    quiver(ax, Xa, Ya, FU(Xa, Ya), FV(Xa, Ya), 1.0, 'Color', 'w', 'LineWidth', 1);
    S.axes(ax); axis(ax, 'equal'); xlim(ax, [0 Q.Lx]); ylim(ax, [0 Q.Ly]); ax.Layer = 'top'; grid(ax, 'off');
    xticks(ax, 0:0.25:1); yticks(ax, 0:0.25:0.75); xlabel(ax, 'Position (x)');
    title(ax, 'Velocity');
    m = F(:,1) <= Q.t(k) + 1e-12;
    ax = nexttile(tl, 5, [1 2]); hold(ax, 'on');
    plot(ax, F(m,4), F(m,5), '-', 'Color', S.blue, 'LineWidth', 2.5, 'DisplayName', 'Finch, 128 x 96');
    plot(ax, E(:,1), E(:,2), 'o', 'Color', S.orange, 'MarkerFaceColor', S.orange, 'MarkerSize', 6, ...
         'DisplayName', 'Martin and Moyce (1952)');
    S.axes(ax); xlim(ax, [0 max(F(:,4))]); ylim(ax, [1 4]);
    xlabel(ax, 'T = t (2g/a)^{1/2}'); ylabel(ax, 'Front position x / a');
    legend(ax, 'Location', 'northwest', 'Box', 'off', 'TextColor', S.ink);
    ax = nexttile(tl, 7, [1 2]); hold(ax, 'on');
    plot(ax, F(m,4), F(m,6), '-', 'Color', S.blue, 'LineWidth', 2.5);
    S.axes(ax); xlim(ax, [0 max(F(:,4))]); ylim(ax, [0 1.05]);
    xlabel(ax, 'T = t (2g/a)^{1/2}'); ylabel(ax, 'Column height h / 2a');
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
