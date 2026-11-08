function s = sweep_series(res, nel, form, it, cfls)
% Max error and status of one form at every CFL number (NaN / "missing" for
% runs not in the table).
    s.err    = NaN(size(cfls));
    s.status = repmat("missing", size(cfls));
    for k = 1:numel(cfls)
        r = res(res.nel == nel & res.form == form & res.iterations == it & ...
                abs(res.CFL - cfls(k)) < 1e-12, :);
        if height(r) == 1
            s.err(k) = r.maxerr;  s.status(k) = r.status;
        end
    end
end
