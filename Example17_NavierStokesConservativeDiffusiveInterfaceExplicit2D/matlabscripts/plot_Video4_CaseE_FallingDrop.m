% Video 4 (case E): the falling drop. Phase field and velocity in time, with
% the drop centroid and speed against free fall drawn up to the frame.
clc; close all; clear;
setenv('CDI_CASE', 'E');
setenv('CDI_FS', '14');
S = plot_style();
Q = load_results(S);
g = 1; y0 = 1.5; R = 0.1; pool = 0.5;
% t, y_c, v_c of the water above the pool (y > pool + 4 eps), from the snapshots
eps = 1 / 96; above = Q.y > pool + 4 * eps;
w = Q.A(above, :); D = [Q.t(:), (w' * Q.y(above)) ./ max(sum(w, 1)', 1e-300), sum(w .* Q.UY(above, :), 1)' ./ max(sum(w, 1)', 1e-300)];
t_imp = sqrt(2 * (y0 - R - pool) / g);
try
    vid = VideoWriter(S.file('Video4_CaseE_FallingDrop.mp4'), 'MPEG-4');
catch
    vid = VideoWriter(S.file('Video4_CaseE_FallingDrop.avi'), 'Motion JPEG AVI');
end
vid.FrameRate = 10; open(vid);
[Xa, Ya] = meshgrid(linspace(0.05, 0.95, 10), linspace(0.05, 1.95, 20));
% CDI_TMAX cuts the movie (a run that blew up after impact); the colour scale
% is set by the frames shown
tmax = str2double(getenv('CDI_TMAX')); if isnan(tmax), tmax = Inf; end
shown = Q.t <= tmax + 1e-9;
umax = max(max(hypot(Q.UX(:, shown), Q.UY(:, shown))));
fig = S.figure(1600, 740);
% CDI_FRAMES limits the number of frames, to check the layout quickly.
nk = nnz(shown);
if ~isempty(getenv('CDI_FRAMES')), nk = min(nk, str2double(getenv('CDI_FRAMES'))); end
for k = 1:nk
    clf(fig);
    tl = tiledlayout(fig, 2, 10, 'TileSpacing', 'compact', 'Padding', 'compact');
    A1 = Q.field(Q.A(:,k));
    ax = nexttile(tl, 1, [2 3]); hold(ax, 'on');
    S.field(ax, Q.Xq, Q.Yq, A1); colormap(ax, S.cmap); clim(ax, [0 1]); add_colorbar(ax, S, 'Phase field (\phi)');
    contour(ax, Q.Xq, Q.Yq, A1, [0.5 0.5], '-', 'LineColor', S.ink, 'LineWidth', 1.5);
    style_field_axes(ax, S, 1, 2); xticks(ax, 0:0.5:1); yticks(ax, 0:0.5:2);
    title(ax, sprintf('Phase field, t = %.2f', Q.t(k)));
    ax = nexttile(tl, 4, [2 3]); hold(ax, 'on');
    S.field(ax, Q.Xq, Q.Yq, Q.field(hypot(Q.UX(:,k), Q.UY(:,k))));
    colormap(ax, S.cmap); clim(ax, [0 umax]); add_colorbar(ax, S, 'Velocity magnitude |u|');
    contour(ax, Q.Xq, Q.Yq, A1, [0.5 0.5], '-', 'LineColor', 'w', 'LineWidth', 1.5);
    FU = scatteredInterpolant(Q.x, Q.y, Q.UX(:,k)); FV = scatteredInterpolant(Q.x, Q.y, Q.UY(:,k));
    quiver(ax, Xa, Ya, FU(Xa, Ya), FV(Xa, Ya), 1.0, 'Color', 'w', 'LineWidth', 1);
    style_field_axes(ax, S, 1, 2); xticks(ax, 0:0.5:1); yticks(ax, 0:0.5:2);
    title(ax, 'Velocity');
    m = D(:,1) <= min(Q.t(k), t_imp) + 1e-12;
    ax = nexttile(tl, 7, [1 4]); hold(ax, 'on');
    plot(ax, D(m,1), D(m,2), '-', 'Color', S.blue, 'LineWidth', 2.5, 'DisplayName', 'Finch, 96 x 192');
    tt = linspace(0, t_imp, 100);
    plot(ax, tt, y0 - g * tt.^2 / 2, '--', 'Color', S.orange, 'LineWidth', 2, 'DisplayName', 'free fall');
    S.axes(ax); xlim(ax, [0 t_imp]); ylim(ax, [pool + R, y0 + 0.05]); ylabel(ax, 'Drop centroid (y)');
    legend(ax, 'Location', 'southwest', 'Box', 'off', 'TextColor', S.ink);
    ax = nexttile(tl, 17, [1 4]); hold(ax, 'on');
    plot(ax, D(m,1), -D(m,3), '-', 'Color', S.blue, 'LineWidth', 2.5);
    plot(ax, tt, g * tt, '--', 'Color', S.orange, 'LineWidth', 2);
    S.axes(ax); xlim(ax, [0 t_imp]); ylim(ax, [0 g * t_imp * 1.05]); xlabel(ax, 'Time (t)'); ylabel(ax, 'Drop speed');
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
