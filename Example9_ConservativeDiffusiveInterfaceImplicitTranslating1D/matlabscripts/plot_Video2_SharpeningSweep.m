% Video 2: one movie per mesh and CFL number from the sharpening-term study,
% written to output/sharpening_sweep/movies/. Columns: 1, 2 and 3 Picard
% iterations. Top: phi near the interface, exact profile dashed. Bottom: the
% error phi - exact. The x axis follows the exact interface, in elements.
% The extrapolated form (one solve per step) is repeated in every column.
clc; close all; clear;
sweep_common;
% Optional: render one mesh or CFL only, e.g. SWEEP_NEL=100 SWEEP_CFL=0.25
if ~isempty(getenv('SWEEP_NEL')); meshes = str2double(getenv('SWEEP_NEL')); end
if ~isempty(getenv('SWEEP_CFL')); cfls = str2double(getenv('SWEEP_CFL')); end
mdir = fullfile(sdir, 'movies');
if ~exist(mdir, 'dir'); mkdir(mdir); end
fps = 8;
win = 15;                                   % half-window, in elements

for nel = meshes
    x   = readmatrix(fullfile(sdir, sprintf('nodes_n%d.csv', nel)));
    dx  = 1/nel;  eps = dx;
    for cfl = cfls
        % every run of this mesh and CFL: column it, the four forms + extrapolated
        runs = struct('it', {}, 'name', {}, 'color', {}, 'ls', {}, 'mk', {}, ...
                      't', {}, 'P', {}, 'status', {});
        for it = iters
            for f = 0:numel(forms)
                if f == 0
                    r = pick(nel, "extrap", 1, cfl); nm = "Extrapolated"; c = xcolor; ls = '--'; mk = 'none';
                else
                    r = pick(nel, forms(f), it, cfl);
                    nm = upper(extractBefore(forms(f), 2)) + extractAfter(forms(f), 1);
                    c = fcolors(f,:); ls = '-'; mk = fmarks(f);
                end
                if height(r) ~= 1; continue; end
                [t, P] = sweep_frames(runfile(r));
                runs(end+1) = struct('it', it, 'name', nm, 'color', c, 'ls', ls, 'mk', mk, ...
                                     't', t, 'P', P, 'status', r.status); %#ok<SAGROW>
            end
        end
        if isempty(runs); continue; end
        nfr = max(arrayfun(@(q) numel(q.t), runs));
        tf  = runs(find(arrayfun(@(q) numel(q.t), runs) == nfr, 1)).t;
        % error range over the whole movie, so the axes do not jump
        emax = 0;
        for q = runs
            for k = 1:numel(q.t)
                xf = x0 + uvel*q.t(k);
                near = abs(x - xf) <= win*dx;
                ex = 0.5*(1 + tanh((x - xf)/(2*eps)));
                emax = max(emax, max(abs(q.P(near,k) - ex(near))));
            end
        end
        emax = min(1.05*emax, 1);

        file = fullfile(mdir, sprintf('Video2_n%d_cfl%g.mp4', nel, cfl));
        vid = VideoWriter(file, 'MPEG-4'); vid.FrameRate = fps; open(vid);
        fig = S.figure(1500, 720);
        for k = 1:nfr
            clf(fig);
            tl = tiledlayout(fig, 2, numel(iters), 'TileSpacing', 'compact', 'Padding', 'compact');
            xf = x0 + uvel*tf(k);
            near = abs(x - xf) <= win*dx;
            xs = (x(near) - xf)/dx;
            ex = 0.5*(1 + tanh((x(near) - xf)/(2*eps)));
            for it = iters
                sel = runs([runs.it] == it);
                axp = nexttile(tl, it); hold(axp, 'on');
                axe = nexttile(tl, numel(iters) + it); hold(axe, 'on');
                plot(axp, xs, ex, '--', 'Color', S.ink, 'LineWidth', 1.2, 'DisplayName', 'Exact');
                h = gobjects(0); lab = strings(0);
                for q = sel
                    kk = min(k, numel(q.t));          % a diverged run holds its last frame
                    pq = q.P(near, kk);
                    eq = pq - 0.5*(1 + tanh((x(near) - x0 - uvel*q.t(kk))/(2*eps)));
                    plot(axp, xs, pq, q.ls, 'Color', q.color, 'LineWidth', 1.4, ...
                         'Marker', q.mk, 'MarkerSize', 4, 'MarkerFaceColor', q.color);
                    h(end+1) = plot(axe, xs, eq, q.ls, 'Color', q.color, 'LineWidth', 1.4, ...
                         'Marker', q.mk, 'MarkerSize', 4, 'MarkerFaceColor', q.color); %#ok<SAGROW>
                    tag = sprintf('%.2g', max(abs(eq)));
                    if q.status == "diverged" && kk < k; tag = 'diverged'; end
                    lab(end+1) = q.name + "  " + tag; %#ok<SAGROW>
                end
                S.axes(axp); S.axes(axe);
                xlim(axp, [-win win]); ylim(axp, [-0.1 1.1]);
                xlim(axe, [-win win]); ylim(axe, [-emax emax]);
                if it == 1
                    title(axp, '1 iteration (lagged = explicit)');
                    ylabel(axp, 'Phase field \phi'); ylabel(axe, 'Error \phi - exact');
                else
                    title(axp, sprintf('%d Picard iterations', it));
                end
                xlabel(axe, '(x - x_{exact}) / \Deltax');
                legend(axe, h, lab, 'Location', 'southwest', 'Box', 'off', ...
                       'TextColor', S.ink2, 'FontSize', S.fs - 2);
            end
            title(tl, sprintf('%d elements, CFL %s, t = %.4f   (legend: current max error)', ...
                  nel, cfllabel(cfl), tf(k)), 'Color', S.ink, 'FontName', S.font, 'FontSize', S.fs + 2);
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
