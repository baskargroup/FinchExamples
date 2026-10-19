% Video 1: animation of the phase field with the numerical interface and the
% exact slotted disk rotated to the current time
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
    % exact outline: rotate the initial outline about (0.5, 0.5) by angle t
    c = cos(times(k)); s = sin(times(k));
    ex = 0.5 + c*(outline(1,:) - 0.5) - s*(outline(2,:) - 0.5);
    ey = 0.5 + s*(outline(1,:) - 0.5) + c*(outline(2,:) - 0.5);
    hE = plot(ax, ex, ey, '--', 'Color', S.orange, 'LineWidth', 2, ...
        'DisplayName', sprintf('Exact (t = %.3f)', times(k)));
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
