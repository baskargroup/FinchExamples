function [U, V] = arrows(Fa, Xa, Ya, ha)
% Velocity u = grad(psi) on the arrow grid, from the interpolant of psi.
    [U, V] = gradient(Fa(Xa, Ya), ha);
end
