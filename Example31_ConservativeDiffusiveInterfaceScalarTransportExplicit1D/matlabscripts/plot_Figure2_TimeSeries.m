clc; close all; clear;
% Figure 2: c and alpha over time, single-hue ramp light (early) -> dark (late),
% laid out as the result panel of Qiu and Mani's 1D demonstration
% alpha_timeseries.csv, c_timeseries.csv: column 1 = x, then one column per snapshot
S     = plot_style();
A     = S.load('alpha_timeseries.csv');
C     = S.load('c_timeseries.csv');
times = S.load('snapshot_times.csv');
x = A(:,1); A = A(:,2:end); C = C(:,2:end);
nsnap = size(A,2);
cmap  = S.ramp(nsnap);

fig = S.figure(800, 600);
tl  = tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
ax1 = nexttile(tl); hold(ax1, 'on');
ax2 = nexttile(tl); hold(ax2, 'on');
for k = 1:nsnap
    plot(ax1, x, C(:,k), '-', 'Color', cmap(k,:), 'LineWidth', 1.5);
    plot(ax2, x, A(:,k), '-', 'Color', cmap(k,:), 'LineWidth', 1.5);
end
for ax = [ax1 ax2]
    S.axes(ax); xlim(ax, [0 1]); xticks(ax, 0:0.1:1);
    colormap(ax, cmap); clim(ax, [times(1) times(end)]);
end
ax1.XTickLabel = [];
ylabel(ax1, 'c / c_{sat}'); ylim(ax1, [-0.1 3.1]);
ylabel(ax2, 'Gas fraction (\alpha)'); ylim(ax2, [-0.05 1.05]);
xlabel(ax2, 'Position (x)');
cb = colorbar(ax2); cb.Layout.Tile = 'east';
cb.Label.String = 'Time'; cb.Label.Color = S.ink;
cb.Color = S.ink2; cb.Box = 'off'; cb.TickDirection = 'out';
S.save(fig, 'Figure2_TimeSeries.pdf');
