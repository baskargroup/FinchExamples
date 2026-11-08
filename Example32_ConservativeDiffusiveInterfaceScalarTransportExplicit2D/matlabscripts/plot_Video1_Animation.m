% Video 1: concentration field with the alpha = 0.5 contour, over time
clc; close all; clear;
plot_common;
try
    vid = VideoWriter(S.file('Video1_Animation.mp4'), 'MPEG-4');          % Mac / Windows
catch
    vid = VideoWriter(S.file('Video1_Animation.avi'), 'Motion JPEG AVI'); % Linux fallback
end
vid.FrameRate = 10;
open(vid);
cmax = max(C(:));
fig = S.figure(620, 520);
for k = 1:nsnap
    clf(fig);
    ax = axes(fig); hold(ax, 'on');
    S.field(ax, Xq, Yq, field(C(:,k)));
    colormap(ax, S.cmap); clim(ax, [0 cmax]);
    cb = colorbar(ax); cb.Label.String = 'c / c_{sat}';
    cb.Label.Color = S.ink; cb.Color = S.ink2; cb.TickDirection = 'in';
    contour(ax, Xq, Yq, field(A(:,k)), [0.5 0.5], '-', 'LineColor', S.ink, 'LineWidth', 1.5);
    style_field_axes(ax, S);
    title(ax, sprintf('t = %.2f', times(k)));
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
