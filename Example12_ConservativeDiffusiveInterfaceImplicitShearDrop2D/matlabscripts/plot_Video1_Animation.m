% Video 1: animation of the phase field with the numerical interface and the
% initial drop, which is the exact solution at t = T
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
    contourf(ax, Xq, Yq, F(Xq, Yq), 20, 'LineColor', 'none');
    colormap(ax, S.ramp(64)); clim(ax, [0 1]);
    cb = colorbar(ax); cb.Label.String = 'Phase field (\phi)';
    cb.Label.Color = S.ink; cb.Color = S.ink2; cb.TickDirection = 'in';
    [~, hN] = contour(ax, Xq, Yq, F(Xq, Yq), [0.5 0.5], '-', ...
        'LineColor', S.ink, 'LineWidth', 1.5, ...
        'DisplayName', sprintf('Numerical (t = %.3f)', times(k)));
    hE = plot(ax, outline(1,:), outline(2,:), '--', 'Color', S.orange, ...
        'LineWidth', 2, 'DisplayName', 'Initial drop (exact at t = T)');
    style_field_axes(ax, S);
    legend(ax, [hE hN], 'Location', 'southwest', 'Box', 'off', 'TextColor', S.ink2);
    drawnow;
    frame = print(fig, '-RGBImage', '-r100');
    if k == 1
        vsz = ceil([size(frame,1) size(frame,2)]/16)*16;
    end
    padded = 255*ones(vsz(1), vsz(2), 3, 'uint8');
    padded(1:size(frame,1), 1:size(frame,2), :) = frame;
    writeVideo(vid, padded);
end
close(vid);
fprintf('video written: %s\n', fullfile(vid.Path, vid.Filename));
