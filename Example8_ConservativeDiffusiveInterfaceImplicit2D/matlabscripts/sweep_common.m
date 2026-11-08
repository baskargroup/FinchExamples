% Shared loading for the sharpening-term study (sharpening_sweep.jl), called
% by plot_Figure6, plot_Figure7 and plot_Video2.
% results.csv:        one row per run (form, iterations, CFL, status, max error ...)
% frames/<run>.bin:   phi at up to 41 instants (read with sweep_frames)
S    = plot_style();
sdir = S.file('sharpening_sweep');
res  = readtable(fullfile(sdir, 'results.csv'), 'TextType', 'string');

% Forms in a fixed order. Color, marker and line style together carry the
% identity, so no series depends on color alone.
forms   = ["lagged" "derived" "reversed" "newton"];
fnames  = ["Lagged \phi^*(1-\phi^*)", "Derived \phi(1-\phi^*)", ...
           "Reversed \phi^*(1-\phi)", "Newton (1-2\phi^*)\phi+\phi^{*2}"];
fcolors = [S.blue; S.orange; [ 22 150 120]/255; [138  79 207]/255];
fmarks  = ['o' 's' '^' 'd'];
xname   = "Extrapolated \phi^e(1-\phi^e), 1 solve";
xcolor  = [107 106 102]/255;

meshes = unique(res.nel)';
cfls   = unique(res.CFL)';
iters  = [1 2 3];

% The row of one run (empty if it is not in the table)
pick = @(nel, form, it, cfl) res(res.nel == nel & res.form == form & ...
                                 res.iterations == it & abs(res.CFL - cfl) < 1e-12, :);
runfile = @(r) fullfile(sdir, 'frames', r.run + ".bin");
cfllabel = @(c) char(replace(string(sprintf('%g', c)), ["0.0625" "0.125" "0.25" "0.5"], ["1/16" "1/8" "1/4" "1/2"]));

% Problem parameters (must match sharpening_sweep.jl)
ux = 1.0;  uy = 0.0;  xc0 = 0.3;  yc0 = 0.5;  R = 0.15;
