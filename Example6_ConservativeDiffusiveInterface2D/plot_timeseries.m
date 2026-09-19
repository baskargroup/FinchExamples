clc;close all; clear;
% Time series of the 2D CDI solution.
% solution_timeseries.csv: columns 1-2 = x, y; columns 3..N+2 = phi at snapshot k
% snapshot_times.csv:      time of each snapshot
ts    = readmatrix('solution_timeseries.csv');
times = readmatrix('snapshot_times.csv');
x     = ts(:,1);
y     = ts(:,2);
PHI   = ts(:,3:end);
nsnap = size(PHI,2);

% Problem parameters (must match main.jl)
ux  = 1.0;  uy = 0.0;
xc0 = 0.3;  yc0 = 0.5;
R   = 0.15;
eps = 1/100;

% Regular grid for plotting
nx = 200; ny = 200;
xq = linspace(min(x), max(x), nx);
yq = linspace(min(y), max(y), ny);
[Xq, Yq] = meshgrid(xq, yq);
F = scatteredInterpolant(x, y, PHI(:,1), 'linear', 'none');
theta = linspace(0, 2*pi, 200);

%% Figure 1: phi = 0.5 interface at every snapshot, colored by time
figure(1);
cmap = parula(nsnap);
hold on;
for k = 1:nsnap
    F.Values = PHI(:,k);
    contour(Xq, Yq, F(Xq, Yq), [0.5 0.5], 'LineColor', cmap(k,:), 'LineWidth', 1.2);
end
xlabel('x'); ylabel('y');
title('2D CDI: interface (\phi = 0.5) over time');
axis image; grid on;
colormap(cmap); cb = colorbar; cb.Label.String = 't';
caxis([times(1) times(end)]);
saveas(gcf, 'cdi2d_timeseries.png');

%% Figure 2: animation (numerical field + exact interface), written to video
fps = 5;
try
    vid = VideoWriter('cdi2d_animation.mp4', 'MPEG-4');          % Mac / Windows
catch
    vid = VideoWriter('cdi2d_animation.avi', 'Motion JPEG AVI'); % Linux fallback
end
vid.FrameRate = fps;
open(vid);

fig2 = figure(2);
set(fig2, 'Position', [100 100 700 600], 'Color', 'w');
for k = 1:nsnap
    clf(fig2);                        % fresh axes every frame (no stacked colorbars)
    F.Values = PHI(:,k);
    xc = xc0 + ux*times(k);
    yc = yc0 + uy*times(k);
    contourf(Xq, Yq, F(Xq, Yq), 20, 'LineColor', 'none'); hold on;
    contour(Xq, Yq, F(Xq, Yq), [0.5 0.5], 'k-', 'LineWidth', 1.5);
    plot(xc + R*cos(theta), yc + R*sin(theta), 'r--', 'LineWidth', 1.5);
    hold off;
    colorbar; caxis([0 1]); axis image;
    xlabel('x'); ylabel('y');
    title(sprintf('2D CDI, t = %.3f  (black: numerical, red: exact)', times(k)));
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

%% Figure 3: total mass (integral of phi) vs time
dx = 1/100;
mass = sum(PHI, 1)' * dx * dx;
figure(3);
plot(times, mass/mass(1), 'bo-');
xlabel('t'); ylabel('mass / initial mass');
title('Mass conservation'); grid on;
saveas(gcf, 'cdi2d_mass.png');
