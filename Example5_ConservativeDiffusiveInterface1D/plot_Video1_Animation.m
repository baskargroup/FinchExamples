% Video 1: animation of numerical vs exact profile, written to a video file
clc; close all; clear;
% solution_timeseries.csv: column 1 = x, columns 2..N+1 = phi at snapshot k
% snapshot_times.csv:      time of each snapshot
S     = plot_style();
ts    = readmatrix('solution_timeseries.csv');
times = readmatrix('snapshot_times.csv');
x     = ts(:,1);
PHI   = ts(:,2:end);
nsnap = size(PHI,2);

% Problem parameters (must match main.jl)
u   = 1.0;
x0  = 0.25;
eps = 1/200;
phi_exact = @(t) 0.5*(1 + tanh((x - x0 - u*t)/(2*eps)));

fps = 5;
try
    vid = VideoWriter('Video1_Animation.mp4', 'MPEG-4');          % Mac / Windows
catch
    vid = VideoWriter('Video1_Animation.avi', 'Motion JPEG AVI'); % Linux fallback
end
vid.FrameRate = fps;
open(vid);

fig2 = S.figure(800, 450);
for k = 1:nsnap
    clf(fig2);
    ax = axes(fig2); hold(ax, 'on');
    % markers first (all nodes, hollow), exact line on top
    hN = plot(ax, x, PHI(:,k), 'o', 'Color', S.blue, 'LineWidth', 1, ...
         'MarkerSize', 5, 'MarkerFaceColor', 'none', ...
         'DisplayName', sprintf('Numerical (t = %.3f)', times(k)));
    hE = plot(ax, x, phi_exact(times(k)), '-', 'Color', S.orange, 'LineWidth', 2, ...
         'DisplayName', sprintf('Exact (t = %.3f)', times(k)));
    S.axes(ax);
    xlabel(ax, 'Position (x)'); ylabel(ax, 'Phase field (\phi)');
    xlim(ax, [0 1]); ylim(ax, [-0.05 1.05]); xticks(ax, 0:0.2:1); yticks(ax, 0:0.25:1);
    grid(ax, 'minor');                       % major + minor grid
    ax.MinorGridColor = S.grid; ax.MinorGridAlpha = 0.6; ax.MinorGridLineStyle = ':';
    ax.XAxis.MinorTickValues = 0:0.05:1;
    ax.YAxis.MinorTickValues = 0:0.05:1;
    legend(ax, [hE hN], 'Location', 'southeast', 'Box', 'off', 'TextColor', S.ink2);   % always empty corner
    drawnow;
    % Render off-screen (getframe can capture a half-drawn window on macOS)
    frame = print(fig2, '-RGBImage', '-r100');
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

