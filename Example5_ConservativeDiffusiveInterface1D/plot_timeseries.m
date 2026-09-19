clc;close all; clear;
% Time series of the 1D CDI solution.
% solution_timeseries.csv: column 1 = x, columns 2..N+1 = phi at snapshot k
% snapshot_times.csv:      time of each snapshot
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

%% Figure 1: all snapshots overlaid
figure(1);
cmap = parula(nsnap);
hold on;
for k = 1:nsnap
    plot(x, PHI(:,k), '-', 'Color', cmap(k,:), 'LineWidth', 1.2);
end
xlabel('x'); ylabel('\phi(x,t)');
title('1D CDI: interface position over time');
ylim([-0.1 1.1]); grid on;
colormap(cmap); cb = colorbar; cb.Label.String = 't';
caxis([times(1) times(end)]);
saveas(gcf, 'cdi1d_timeseries.png');

%% Figure 2: animation (numerical vs exact), also written to a video file
fps = 5;                              % frames per second in the video
try
    vid = VideoWriter('cdi1d_animation.mp4', 'MPEG-4');   % Mac / Windows
catch
    vid = VideoWriter('cdi1d_animation.avi', 'Motion JPEG AVI');  % Linux fallback
end
vid.FrameRate = fps;
open(vid);

fig2 = figure(2);
set(fig2, 'Position', [100 100 800 450], 'Color', 'w');   % fixed size for video
for k = 1:nsnap
    clf(fig2);                        % fresh axes every frame
    plot(x, PHI(:,k), 'bo-', 'MarkerSize', 3, 'DisplayName', 'Numerical'); hold on;
    plot(x, phi_exact(times(k)), 'r-', 'DisplayName', 'Exact'); hold off;
    xlabel('x'); ylabel('\phi');
    title(sprintf('1D CDI, t = %.4f', times(k)));
    ylim([-0.1 1.1]); grid on; legend('Location', 'northwest');
    drawnow;
    % Render off-screen instead of grabbing the screen (getframe can capture
    % a half-drawn window on macOS and produce striped frames)
    frame = print(fig2, '-RGBImage', '-r100');
    % Encoders want every frame the same size, with even (ideally x16)
    % dimensions; crop to a fixed size decided from the first frame.
    if k == 1
        vsz = floor([size(frame,1) size(frame,2)]/16)*16;
    end
    frame = frame(1:vsz(1), 1:vsz(2), :);
    writeVideo(vid, frame);
end
close(vid);
fprintf('video written: %s\n', fullfile(vid.Path, vid.Filename));

%% Figure 3: interface position (phi = 0.5 crossing) vs time
xi = zeros(nsnap,1);
for k = 1:nsnap
    [~, imax] = max(PHI(:,k) >= 0.5);   % first node with phi >= 0.5
    xi(k) = x(imax);
end
figure(3);
plot(times, xi, 'bo', 'DisplayName', 'Numerical'); hold on;
plot(times, x0 + u*times, 'r-', 'DisplayName', 'Exact x_0 + u t');
xlabel('t'); ylabel('interface position');
title('Interface position vs time'); grid on; legend('Location', 'northwest');
saveas(gcf, 'cdi1d_interface_position.png');
