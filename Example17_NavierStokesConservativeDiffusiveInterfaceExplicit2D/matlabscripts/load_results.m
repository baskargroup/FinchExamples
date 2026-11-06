function R = load_results(S, folder)
% Read one run of main.jl. folder: results folder relative to the example
% (default: the case's folder from plot_style). Fields on a regular grid:
%   R.Xq, R.Yq     plotting grid (about 200 points per unit length)
%   R.field(v)     nodal values v -> values on the grid
%   R.A, R.UX, R.UY, R.P   nodal values, one column per snapshot
%   R.t            snapshot times
%   R.h            history.csv columns: t, volume, kinetic energy,
%                  ||div u||, centroid x, y, velocity x, y, circularity,
%                  ||u - u0||
    if nargin < 2, d = S.outdir; else, d = S.dir(folder); end
    rd = @(n) readmatrix(fullfile(d, n));
    At = rd('alpha_timeseries.csv'); Ut = rd('ux_timeseries.csv');
    Vt = rd('uy_timeseries.csv');    Pt = rd('p_timeseries.csv');
    R.x = At(:,1); R.y = At(:,2);
    R.A = At(:,3:end); R.UX = Ut(:,3:end); R.UY = Vt(:,3:end); R.P = Pt(:,3:end);
    R.t = rd('snapshot_times.csv');
    R.h = rd('history.csv');
    R.Lx = max(R.x); R.Ly = max(R.y);
    [R.Xq, R.Yq] = meshgrid(linspace(0, R.Lx, round(200*R.Lx)), ...
                            linspace(0, R.Ly, round(200*R.Ly)));
    F = scatteredInterpolant(R.x, R.y, R.A(:,1), 'linear', 'nearest');
    R.field = @(v) interp_values(F, v, R.Xq, R.Yq);
    R.dir = d;
end

function out = interp_values(F, v, X, Y)
    F.Values = v;
    out = F(X, Y);
end
