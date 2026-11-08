% Video 2: gas fraction field with the alpha = 0.5 contour of the initial
% bubbles (dashed) and the current ones (solid), over time
clc; close all; clear;
plot_common;
try
    vid = VideoWriter(S.file('Video2_GasFraction.mp4'), 'MPEG-4');          % Mac / Windows
catch
    vid = VideoWriter(S.file('Video2_GasFraction.avi'), 'Motion JPEG AVI'); % Linux fallback
end
vid.FrameRate = 10;
open(vid);
A0 = field(A(:,1));
fig = S.figure(620, 520);
for k = 1:nsnap
    clf(fig);
    ax = axes(fig); hold(ax, 'on');
    Ak = field(A(:,k));
    S.field(ax, Xq, Yq, Ak);
    colormap(ax, S.cmap); clim(ax, [0 1]);
    cb = colorbar(ax); cb.Label.String = 'Gas fraction (\alpha)';
    cb.Label.Color = S.ink; cb.Color = S.ink2; cb.TickDirection = 'in';
    [~, h0] = contour(ax, Xq, Yq, A0, [0.5 0.5], '--', 'LineColor', 'w', ...
                      'LineWidth', 1.5, 'DisplayName', 'Initial');
    [~, hk] = contour(ax, Xq, Yq, Ak, [0.5 0.5], '-', 'LineColor', S.ink, ...
                      'LineWidth', 1.5, 'DisplayName', 'Current');
    style_field_axes(ax, S);
    legend(ax, [h0 hk], 'Location', 'southeast', 'Box', 'off', 'TextColor', 'w');   % on the dark-blue liquid
    title(ax, sprintf('t = %.2f', times(k)));
    drawnow;
    % Render off-screen (getframe can capture a half-drawn window on macOS)
    frame = print(fig, '-RGBImage', '-r100');
    if k == 1
        vsz = ceil([size(frame,1) size(frame,2)]/16)*16;    % encoder-friendly size
    end
    % pad with white (never crop, so no labels are lost) to a fixed size
    padded = 255*ones(vsz(1), vsz(2), 3, 'uint8');
    padded(1:size(frame,1), 1:size(frame,2), :) = frame;
    writeVideo(vid, padded);
end
close(vid);
fprintf('video written: %s\n', fullfile(vid.Path, vid.Filename));
