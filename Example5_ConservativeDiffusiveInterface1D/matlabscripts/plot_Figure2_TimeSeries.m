% Figure 2: interface profiles over time, single-hue ramp light (early) -> dark (late)
clc; close all; clear;
% solution_timeseries.csv: column 1 = x, columns 2..N+1 = phi at snapshot k
% snapshot_times.csv:      time of each snapshot
S     = plot_style();
ts    = S.load('solution_timeseries.csv');
times = S.load('snapshot_times.csv');
x     = ts(:,1);
PHI   = ts(:,2:end);
nsnap = size(PHI,2);

% Problem parameters (must match main.jl)
u   = 1.0;
x0  = 0.25;
eps = 1/200;

cmap = S.ramp(nsnap);
fig1 = S.figure(800, 450);
ax = axes(fig1); hold(ax, 'on');
for k = 1:nsnap
    plot(ax, x, PHI(:,k), '-', 'Color', cmap(k,:), 'LineWidth', 1.5);
end
S.axes(ax);
xlabel(ax, 'Position (x)'); ylabel(ax, 'Phase field (\phi)');
xlim(ax, [0.15 0.6]); ylim(ax, [-0.05 1.05]); yticks(ax, 0:0.25:1);
grid(ax, 'minor');                       % major + minor grid
ax.MinorGridColor = S.grid; ax.MinorGridAlpha = 0.6; ax.MinorGridLineStyle = ':';
ax.XAxis.MinorTickValues = 0:0.01:1;
ax.YAxis.MinorTickValues = 0:0.05:1;
colormap(ax, cmap); clim(ax, [times(1) times(end)]);
cb = colorbar(ax); cb.Label.String = 'Time'; cb.Label.Color = S.ink;
cb.Color = S.ink2; cb.Box = 'off'; cb.TickDirection = 'out';
S.save(fig1, 'Figure2_TimeSeries.pdf');

