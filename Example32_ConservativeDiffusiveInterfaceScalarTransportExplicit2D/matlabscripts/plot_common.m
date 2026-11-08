% Shared loading for the plot scripts (called by each script).
% alpha_timeseries.csv, c_timeseries.csv (and psi_timeseries.csv for the
% electrode): columns 1-2 = x, y; then one column per snapshot.
% snapshot_times.csv: time of each snapshot. The case comes from a variable
% casename set by the calling script, else CDI_CASE, else 'bubbles'.
if ~exist('casename', 'var'), casename = ''; end
S     = plot_style(casename);
ELECTRODE = strcmp(S.casename, 'electrode');
At    = S.load('alpha_timeseries.csv');
Ct    = S.load('c_timeseries.csv');
times = S.load('snapshot_times.csv');
x = At(:,1); y = At(:,2);
A = At(:,3:end); C = Ct(:,3:end);
nsnap = size(A,2);

% Regular grid for contour plotting (the periodic electrode box extrapolates
% to the nearest node at the sides)
nq = 400;
[Xq, Yq] = meshgrid(linspace(0, 1, nq), linspace(0, 1, nq));
if ELECTRODE
    F = scatteredInterpolant(x, y, A(:,1), 'linear', 'nearest');
else
    F = scatteredInterpolant(x, y, A(:,1), 'linear', 'none');
end
field = @(v) subsref(setfield(F, 'Values', v), struct('type', '()', 'subs', {{Xq, Yq}}));

if ELECTRODE
    Pt = S.load('psi_timeseries.csv');
    P  = Pt(:,3:end);
    % Coarse grid for the velocity arrows, u = grad(psi)
    na = 13; ha = 1/(na - 1);
    [Xa, Ya] = meshgrid(linspace(0, 1, na), linspace(0, 1, na));
    Fa = scatteredInterpolant(x, y, P(:,1), 'linear', 'nearest');
    velocity = @(v) arrows(setfield(Fa, 'Values', v), Xa, Ya, ha);
    % c colour range: saturated (1) up to the largest value reached
    cmax = max(C(:));
end
