function S = plot_style()
% Shared figure styling for the Finch examples.
% Returns a struct with colors and helper functions:
%   S.blue, S.orange       series colors (numerical, exact)
%   S.ink, S.ink2          text colors (primary, secondary)
%   S.grid                 grid line color
%   S.ramp(n)              n-step single-hue blue ramp (light -> dark) for time
%   S.cmap                 colormap of the 2D fields (jet, as in the slides)
%   S.field(ax, X, Y, F)   draw a field on a grid, colours interpolated per cell
%   S.figure(w, h)         new white figure of w x h pixels
%   S.axes(ax)             apply axis styling
%   S.load(name)           read a result file from the case's output folder
%   S.file(name)           path of a file in that folder
%   S.save(fig, name)      export as PDF into that folder (vector, or 300 dpi with fields)
S.blue   = [ 42 120 214]/255;   % #2a78d6
S.orange = [235 104  52]/255;   % #eb6834
S.ink    = [ 11  11  11]/255;   % #0b0b0b
S.ink2   = [ 82  81  78]/255;   % #52514e
S.grid   = [220 219 214]/255;
S.font   = 'Helvetica';
S.fs     = 12;
% CDI_FS overrides the font size, for figures shown large on a slide.
fs = getenv('CDI_FS'); if ~isempty(fs), S.fs = str2double(fs); end

ramp_ends = [205 226 251; 13 54 107]/255;   % #cde2fb -> #0d366b
S.ramp = @(n) [linspace(ramp_ends(1,1), ramp_ends(2,1), n)', ...
               linspace(ramp_ends(1,2), ramp_ends(2,2), n)', ...
               linspace(ramp_ends(1,3), ramp_ends(2,3), n)'];

% 2D fields in the look of the slides: jet, smooth (no contour bands). The
% surface sits at z = -1 so contour lines and arrows drawn after it stay on top.
S.cmap  = jet(256);
S.field = @(ax, X, Y, F) surface(ax, X, Y, -ones(size(F)), F, ...
                                 'EdgeColor', 'none', 'FaceColor', 'interp');

S.figure = @(w, h) figure('Color', 'w', 'Position', [100 100 w h], ...
                          'InvertHardcopy', 'off');
S.axes   = @style_axes;
% Results live in the example's output_case<A|B|C>/ folder, one level up from
% this script. The case comes from CDI_CASE (default A); CDI_OUTDIR (relative
% to the example folder) selects any other results folder.
here = fileparts(mfilename('fullpath'));
cs = getenv('CDI_CASE'); if isempty(cs), cs = 'A'; end
outdir = fullfile(here, '..', ['output_case' upper(cs)]);
if ~isempty(getenv('CDI_OUTDIR'))
    outdir = fullfile(here, '..', getenv('CDI_OUTDIR'));
end
S.outdir = outdir;
S.dir    = @(name) fullfile(here, '..', name);   % another results folder
S.file   = @(name) fullfile(outdir, name);
S.load   = @(name) readmatrix(fullfile(outdir, name));
% Vector PDF, except figures with a smooth field (S.field): a shaded surface
% in a vector PDF shows seams, so those are written as a 300 dpi image.
S.save   = @(fig, name) save_figure(fig, fullfile(outdir, name));

    function style_axes(ax)
        set(ax, 'FontName', S.font, 'FontSize', S.fs, ...
                'XColor', S.ink2, 'YColor', S.ink2, ...
                'LineWidth', 0.8, 'Box', 'on', 'TickDir', 'in', ...   % full box, ticks inside
                'TickLength', [0.01 0.01], ...
                'GridColor', S.grid, 'GridAlpha', 1, 'Layer', 'bottom');
        grid(ax, 'on');
        ax.Toolbar.Visible = 'off';     % keep the toolbar out of exported images
        ax.Title.Color = S.ink;  ax.Title.FontWeight = 'normal';
        ax.Title.FontSize = S.fs + 2;
        ax.XLabel.Color = S.ink; ax.YLabel.Color = S.ink;
    end
end

function save_figure(fig, file)
    % a 200 dpi PNG next to every PDF, for the slides
    [d, n] = fileparts(file);
    exportgraphics(fig, fullfile(d, [n '.png']), 'Resolution', 200, 'BackgroundColor', 'w');
    if isempty(findobj(fig, 'Type', 'surface'))
        exportgraphics(fig, file, 'ContentType', 'vector', 'BackgroundColor', 'w');
    else
        exportgraphics(fig, file, 'ContentType', 'image', 'Resolution', 300, ...
                       'BackgroundColor', 'w');
    end
end
