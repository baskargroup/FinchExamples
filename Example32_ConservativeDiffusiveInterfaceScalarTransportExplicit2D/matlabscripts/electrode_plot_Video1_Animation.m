% Video 1: gas fraction | dissolved gas with the bubble outline and the flow
clc; close all; clear;
casename = 'electrode';
plot_common;
try
    vid = VideoWriter(S.file('Video1_Animation.mp4'), 'MPEG-4');          % Mac / Windows
catch
    vid = VideoWriter(S.file('Video1_Animation.avi'), 'Motion JPEG AVI'); % Linux fallback
end
vid.FrameRate = 12;
open(vid);
fig = S.figure(1200, 520);
for k = 1:nsnap
    clf(fig);
    tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    [U, V] = velocity(P(:,k));
    draw_frame(tl, S, Xq, Yq, field(A(:,k)), field(C(:,k)), Xa, Ya, U, V, cmax, times(k));
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
