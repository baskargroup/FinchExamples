% Video 2: one movie per mesh and CFL number from the sharpening-term study,
% written to output/sharpening_sweep/movies/. Each panel is the error
% phi - exact in a window that follows the exact circle, with the exact
% interface (dashed) and the computed phi = 0.5 contour (solid). Rows: 1, 2
% and 3 Picard iterations; columns: the forms. The extrapolated form takes one
% solve per step and appears in the first row only. The color scale is the
% same in every movie and saturates at +-elim.
clc; close all; clear;
sweep_common;
% Optional: render one mesh or CFL only, e.g. SWEEP_NEL=100 SWEEP_CFL=0.25
if ~isempty(getenv('SWEEP_NEL')); meshes = str2double(getenv('SWEEP_NEL')); end
if ~isempty(getenv('SWEEP_CFL')); cfls = str2double(getenv('SWEEP_CFL')); end
mdir = fullfile(sdir, 'movies');
if ~exist(mdir, 'dir'); mkdir(mdir); end
fps  = 8;
half = 0.25;                                % window half-width around the centre
elim = 0.03;                                % color limit of the error
cols = ["extrap" forms];
ctitle = ["Extrapolated" "Lagged" "Derived" "Reversed" "Newton"];

% diverging map: blue (phi too low) - light neutral - orange (phi too high)
mid  = [0.95 0.95 0.94];
cmap = [interp1([0 1], [S.blue; mid], linspace(0, 1, 128)); ...
        interp1([0 1], [mid; S.orange], linspace(0, 1, 128))];
theta = linspace(0, 2*pi, 400);

for nel = meshes
    xy  = readmatrix(fullfile(sdir, sprintf('nodes_n%d.csv', nel)));
    dx  = 1/nel;  eps = dx;
    % nodes of the quad mesh on the (nel+1)^2 grid: G(j,i) is the node at
    % x = (i-1)*dx, y = (j-1)*dx
    ix = round(xy(:,1)/dx) + 1;  iy = round(xy(:,2)/dx) + 1;
    G  = zeros(nel+1);  G(sub2ind(size(G), iy, ix)) = 1:size(xy,1);
    gv = (0:nel)*dx;
    for cfl = cfls
        runs = struct('row', {}, 'col', {}, 't', {}, 'P', {}, 'status', {}, 'tend', {});
        for c = 1:numel(cols)
            for it = iters
                if cols(c) == "extrap" && it ~= 1; continue; end
                r = pick(nel, cols(c), it, cfl);
                if height(r) ~= 1; continue; end
                [t, P] = sweep_frames(runfile(r));
                runs(end+1) = struct('row', it, 'col', c, 't', t, 'P', P, 'status', r.status, 'tend', r.t_end); %#ok<SAGROW>
            end
        end
        if isempty(runs); continue; end
        nfr = max(arrayfun(@(q) numel(q.t), runs));
        tf  = runs(find(arrayfun(@(q) numel(q.t), runs) == nfr, 1)).t;

        file = fullfile(mdir, sprintf('Video2_n%d_cfl%g.mp4', nel, cfl));
        vid = VideoWriter(file, 'MPEG-4'); vid.FrameRate = fps; open(vid);
        fig = S.figure(1500, 880);
        for k = 1:nfr
            clf(fig);
            tl = tiledlayout(fig, numel(iters), numel(cols), 'TileSpacing', 'compact', 'Padding', 'compact');
            for q = runs
                ax = nexttile(tl, (q.row - 1)*numel(cols) + q.col); hold(ax, 'on');
                kk = min(k, numel(q.t));            % a diverged run holds its last frame
                tq = q.t(kk);
                xc = xc0 + ux*tq;  yc = yc0 + uy*tq;
                isx = find(abs(gv - xc) <= half);  isy = find(abs(gv - yc) <= half);
                [X, Y] = meshgrid(gv(isx), gv(isy));
                PH = q.P(G(isy, isx), kk);  PH = reshape(PH, size(X));
                EX = 0.5*(1 - tanh((sqrt((X - xc).^2 + (Y - yc).^2) - R)/(2*eps)));
                imagesc(ax, gv(isx), gv(isy), PH - EX, [-elim elim]);
                colormap(ax, cmap);
                plot(ax, xc + R*cos(theta), yc + R*sin(theta), '--', 'Color', S.ink, 'LineWidth', 1.0);
                contour(ax, X, Y, PH, [0.5 0.5], '-', 'LineColor', S.ink, 'LineWidth', 1.0);
                set(ax, 'YDir', 'normal', 'XTick', [], 'YTick', [], 'Box', 'on', 'Layer', 'top', ...
                        'XColor', S.ink2, 'YColor', S.ink2);
                axis(ax, 'equal');
                xlim(ax, xc + [-half half]); ylim(ax, yc + [-half half]);
                emax = max(abs(PH(:) - EX(:)));
                if q.status == "diverged" && kk < k
                    tag = sprintf('diverged at t = %.3f', q.tend);
                else
                    tag = sprintf('max error %.2g', emax);
                end
                if q.col == 1
                    head = 'Extrapolated, 1 solve';
                elseif q.col == 2 && q.row == 1
                    head = 'Lagged, 1 it. (explicit)';
                else
                    head = sprintf('%s, %d it.', ctitle(q.col), q.row);
                end
                title(ax, {head, tag}, 'FontName', S.font, 'FontSize', S.fs - 1, ...
                      'FontWeight', 'normal', 'Color', S.ink);
                ax.Toolbar.Visible = 'off';
            end
            % the empty tiles under the extrapolated form carry the key
            ax = nexttile(tl, numel(cols) + 1); axis(ax, 'off');
            colormap(ax, cmap); clim(ax, [-elim elim]);
            cb = colorbar(ax, 'west'); cb.Label.String = 'Error \phi - exact';
            cb.Color = S.ink2; cb.Label.Color = S.ink; cb.FontSize = S.fs - 1;
            ax = nexttile(tl, 2*numel(cols) + 1); axis(ax, 'off');
            text(ax, 0, 0.6, {'dashed: exact interface', 'solid: computed \phi = 0.5', ...
                 sprintf('window: centre \\pm %.2f', half), sprintf('colors saturate at \\pm%.2f', elim)}, ...
                 'FontName', S.font, 'FontSize', S.fs - 1, 'Color', S.ink2, 'VerticalAlignment', 'middle');
            title(tl, sprintf('%d^2 elements, CFL %s, t = %.4f', nel, cfllabel(cfl), tf(k)), ...
                  'Color', S.ink, 'FontName', S.font, 'FontSize', S.fs + 2);
            drawnow;
            frame = print(fig, '-RGBImage', '-r100');
            if k == 1; vsz = ceil([size(frame,1) size(frame,2)]/16)*16; end
            padded = 255*ones(vsz(1), vsz(2), 3, 'uint8');
            padded(1:min(end,size(frame,1)), 1:min(end,size(frame,2)), :) = ...
                frame(1:min(vsz(1),size(frame,1)), 1:min(vsz(2),size(frame,2)), :);
            writeVideo(vid, padded);
        end
        close(vid); close(fig);
        fprintf('video written: %s\n', file);
    end
end
