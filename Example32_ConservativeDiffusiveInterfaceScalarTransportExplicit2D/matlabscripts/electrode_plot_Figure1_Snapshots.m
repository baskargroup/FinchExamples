% Figure 1: gas fraction and dissolved gas (with the flow) at three times
clc; close all; clear;
casename = 'electrode';
plot_common;
fig = S.figure(1100, 1350);
tl  = tiledlayout(fig, 3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for f = [1/3 2/3 1]
    [~, k] = min(abs(times - f*times(end)));
    [U, V] = velocity(P(:,k));
    draw_frame(tl, S, Xq, Yq, field(A(:,k)), field(C(:,k)), Xa, Ya, U, V, cmax, times(k));
end
S.save(fig, 'Figure1_Snapshots.pdf');
