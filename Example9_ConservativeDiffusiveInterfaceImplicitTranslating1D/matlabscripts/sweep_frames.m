function [t, P] = sweep_frames(file)
% Read the frames written by sharpening_sweep.jl: Int64 nnode, Int64 nf,
% Float64 times(nf), Float32 phi(nnode, nf). A diverged run has fewer frames.
    fid = fopen(file, 'r');
    n  = fread(fid, 1, 'int64');
    nf = fread(fid, 1, 'int64');
    t  = fread(fid, nf, 'float64');
    P  = fread(fid, [n nf], 'float32=>double');
    fclose(fid);
end
