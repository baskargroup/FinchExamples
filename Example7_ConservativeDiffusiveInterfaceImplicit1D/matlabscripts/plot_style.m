function S = plot_style()
% Shared figure styling for the Finch examples.
% Returns a struct with colors and helper functions:
%   S.blue, S.orange       series colors (numerical, exact)
%   S.ink, S.ink2          text colors (primary, secondary)
%   S.grid                 grid line color
%   S.ramp(n)              n-step single-hue blue ramp (light -> dark) for time
%   S.figure(w, h)         new white figure of w x h pixels
%   S.axes(ax)             apply axis styling
%   S.load(name)           read a result file from output/
%   S.file(name)           path of a file in output/
%   S.save(fig, name)      export as vector PDF into output/
S.blue   = [ 42 120 214]/255;   % #2a78d6
S.orange = [235 104  52]/255;   % #eb6834
S.ink    = [ 11  11  11]/255;   % #0b0b0b
S.ink2   = [ 82  81  78]/255;   % #52514e
S.grid   = [220 219 214]/255;
S.font   = 'Helvetica';
S.fs     = 12;

ramp_ends = [205 226 251; 13 54 107]/255;   % #cde2fb -> #0d366b
S.ramp = @(n) [linspace(ramp_ends(1,1), ramp_ends(2,1), n)', ...
               linspace(ramp_ends(1,2), ramp_ends(2,2), n)', ...
               linspace(ramp_ends(1,3), ramp_ends(2,3), n)'];

S.figure = @(w, h) figure('Color', 'w', 'Position', [100 100 w h], ...
                          'InvertHardcopy', 'off');
S.axes   = @style_axes;
% Results live in the example's output/ folder, one level up from this
% script. Deriving it from the script's own path keeps these helpers
% working whatever the current directory happens to be.
outdir   = fullfile(fileparts(mfilename('fullpath')), '..', 'output');
S.file   = @(name) fullfile(outdir, name);
S.load   = @(name) readmatrix(fullfile(outdir, name));
S.save   = @(fig, name) exportgraphics(fig, fullfile(outdir, name), ...
                                       'ContentType', 'vector', ...
                                       'BackgroundColor', 'w');

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
