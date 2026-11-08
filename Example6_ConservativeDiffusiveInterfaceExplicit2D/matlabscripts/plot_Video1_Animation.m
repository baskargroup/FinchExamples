% Video 1: animation of the phase field with numerical and exact interfaces
clc; close all; clear;
plot_common;
fps = 5;
try
    vid = VideoWriter(S.file('Video1_Animation.mp4'), 'MPEG-4');          % Mac / Windows
catch
    vid = VideoWriter(S.file('Video1_Animation.avi'), 'Motion JPEG AVI'); % Linux fallback
end
vid.FrameRate = fps;
open(vid);

fig = S.figure(620, 520);
for k = 1:nsnap
    clf(fig);
    ax = axes(fig); hold(ax, 'on');
    F.Values = PHI(:,k);
    xc = xc0 + ux*times(k);  yc = yc0 + uy*times(k);
    contourf(ax, Xq, Yq, F(Xq, Yq), 20, 'LineColor', 'none');
    colormap(ax, S.ramp(64)); clim(ax, [0 1]);
    cb = colorbar(ax); cb.Label.String = 'Phase field (\phi)';
    cb.Label.Color = S.ink; cb.Color = S.ink2; cb.TickDirection = 'in';
    [~, hN] = contour(ax, Xq, Yq, F(Xq, Yq), [0.5 0.5], '-', ...
        'LineColor', S.ink, 'LineWidth', 1.5, ...
        'DisplayName', sprintf('Numerical (t = %.3f)', times(k)));
    hE = plot(ax, xc + R*cos(theta), yc + R*sin(theta), '--', ...
        'Color', S.orange, 'LineWidth', 2, ...
        'DisplayName', sprintf('Exact (t = %.3f)', times(k)));
    style_field_axes(ax, S);
    legend(ax, [hE hN], 'Location', 'northwest', 'Box', 'off', 'TextColor', S.ink2);
    drawnow;
    % Render off-screen (getframe can capture a half-drawn window on macOS)
    frame = print(fig, '-RGBImage', '-r100');
    if k == 1
        vsz = ceil([size(frame,1) size(frame,2)]/16)*16;    % encoder-friendly size
    end
    padded = 255*ones(vsz(1), vsz(2), 3, 'uint8');          % pad, never crop
    padded(1:size(frame,1), 1:size(frame,2), :) = frame;
    writeVideo(vid, padded);
end
close(vid);
fprintf('video written: %s\n', fullfile(vid.Path, vid.Filename));
