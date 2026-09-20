%%
clc;clear;close all

% Load data
outdir = fullfile(fileparts(mfilename('fullpath')), '..', 'output');
data = readmatrix(fullfile(outdir, 'solution.csv'));
x = data(:,1);
y = data(:,2);
u = data(:,3);

% Interpolate onto a regular grid for plotting
nx = 50;
ny = 50;
xq = linspace(min(x), max(x), nx);
yq = linspace(min(y), max(y), ny);
[Xq, Yq] = meshgrid(xq, yq);

F = scatteredInterpolant(x, y, u, 'linear', 'none');
Uq = F(Xq, Yq);

% 2D contour plot - numerical
figure(1);
contourf(Xq, Yq, Uq, 20, 'LineColor', 'none');
colorbar;
xlabel('x');
ylabel('y');
title('2D Convection-Diffusion Solution (Numerical)');
axis image;

% Exact solution contour
Uexact = sin(pi*Xq).*sin(pi*Yq);
figure(2);
contourf(Xq, Yq, Uexact, 20, 'LineColor', 'none');
colorbar;
xlabel('x');
ylabel('y');
title('2D Convection-Diffusion Solution (Exact)');
axis image;

% Pointwise error
figure(3);
contourf(Xq, Yq, abs(Uq - Uexact), 20, 'LineColor', 'none');
colorbar;
xlabel('x');
ylabel('y');
title('Absolute Error');
axis image;

fprintf('Max nodal error (interpolated grid): %.4e\n', max(abs(Uq(:) - Uexact(:))));