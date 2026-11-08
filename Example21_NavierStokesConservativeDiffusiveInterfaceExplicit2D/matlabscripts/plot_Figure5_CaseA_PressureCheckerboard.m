% Figure 5 (case A): pressure of the drop advection run with the consistent
% term (output_caseA) at the end. The velocity error of that run grows in
% time; this shows the pressure carries a node-to-node (checkerboard)
% oscillation of the same size as its smooth part, that the velocity error
% has the same pattern, and that both grow at the same rate.
clc; close all; clear;
setenv('CDI_CASE', 'A');
setenv('CDI_FS', '14');
S = plot_style(); Q = load_results(S, 'output_caseA');
% nodal grid: 64 x 64
[xs,~,ix] = unique(Q.x); [ys,~,iy] = unique(Q.y); nx = numel(xs); ny = numel(ys);
G = @(v) accumarray([iy ix], v, [ny nx]);
k = numel(Q.t);  P = G(Q.P(:,k));  A = G(Q.A(:,k));
E = G(hypot(Q.UX(:,k)-1, Q.UY(:,k)));
% checkerboard measure: p minus mean of 4 periodic neighbours
Pn = (circshift(P,1,1)+circshift(P,-1,1)+circshift(P,1,2)+circshift(P,-1,2))/4;
C = P - Pn;
fprintf('t = %.3f  p range [%g, %g]  max|p - nbr avg| = %g  max vel err = %g\n', Q.t(k), min(P(:)), max(P(:)), max(abs(C(:))), max(E(:)));
% time history of checkerboard measure and velocity error
nt = numel(Q.t); cb = zeros(nt,1); ve = zeros(nt,1); pr = zeros(nt,1);
for j = 1:nt
  Pj = G(Q.P(:,j)); Pnj = (circshift(Pj,1,1)+circshift(Pj,-1,1)+circshift(Pj,1,2)+circshift(Pj,-1,2))/4;
  cb(j) = max(abs(Pj(:)-Pnj(:))); pr(j) = max(Pj(:))-min(Pj(:));
  ve(j) = max(hypot(Q.UX(:,j)-1, Q.UY(:,j)));
end
fig = S.figure(1400, 1000); tl = tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');
ax = nexttile(tl); S.field(ax, xs, ys, P); colormap(ax, S.cmap); colorbar(ax); hold(ax,'on');
contour(ax, xs, ys, A, [0.5 0.5], 'w', 'LineWidth', 1.5); clim(ax, max(abs(P(:)))*[-1 1]); axis(ax,'equal','tight'); S.axes(ax); title(ax, sprintf('p at t = %.2f', Q.t(k)));
ax = nexttile(tl); S.field(ax, xs, ys, C); colormap(ax, S.cmap); colorbar(ax); hold(ax,'on');
contour(ax, xs, ys, A, [0.5 0.5], 'w', 'LineWidth', 1.5); clim(ax, max(abs(C(:)))*[-1 1]); axis(ax,'equal','tight'); S.axes(ax); title(ax, 'p - mean of 4 neighbours');
ax = nexttile(tl); S.field(ax, xs, ys, E); colormap(ax, S.cmap); colorbar(ax); hold(ax,'on');
contour(ax, xs, ys, A, [0.5 0.5], 'w', 'LineWidth', 1.5); clim(ax, [0 max(E(:))]); axis(ax,'equal','tight'); S.axes(ax); title(ax, '|u - (1,0)|');
jm = find(abs(ys-0.5) == min(abs(ys-0.5)), 1);
ax = nexttile(tl); plot(ax, xs, P(jm,:), '.-', 'Color', S.blue, 'LineWidth', 1.2); S.axes(ax); title(ax, sprintf('p along y = %.3f', ys(jm))); xlabel(ax,'x');
ax = nexttile(tl); plot(ax, xs, E(jm,:), '.-', 'Color', S.orange, 'LineWidth', 1.2); S.axes(ax); title(ax, '|u - (1,0)| along the same line'); xlabel(ax,'x');
ax = nexttile(tl); semilogy(ax, Q.t, max(cb,1e-16), '-', 'Color', S.blue, 'LineWidth', 2, 'DisplayName', 'max |p - nbr avg|'); hold(ax,'on');
semilogy(ax, Q.t, max(pr,1e-16), '--', 'Color', S.ink2, 'LineWidth', 1.5, 'DisplayName', 'max p - min p');
semilogy(ax, Q.t, max(ve,1e-16), '-', 'Color', S.orange, 'LineWidth', 2, 'DisplayName', 'max |u - (1,0)|');
S.axes(ax); legend(ax, 'Location', 'southeast', 'Box', 'off'); xlabel(ax, 't'); title(ax, 'in time');
S.save(fig, 'Figure5_CaseA_PressureCheckerboard.pdf');
fprintf('written: %s\n', S.file('Figure5_CaseA_PressureCheckerboard.pdf'));
