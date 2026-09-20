%% 
clc;clear;close all

% Load data
outdir = fullfile(fileparts(mfilename('fullpath')), '..', 'output');
data = readmatrix(fullfile(outdir, 'solution.csv'));
x = data(:,1);
y = data(:,2);
u = data(:,3);

% Interpolate onto a regular grid for plotting
nx = 50; % number of grid points in x
ny = 50; % number of grid points in y
xq = linspace(min(x), max(x), nx);
yq = linspace(min(y), max(y), ny);
[Xq, Yq] = meshgrid(xq, yq);

F = scatteredInterpolant(x, y, u, 'linear', 'none');
Uq = F(Xq, Yq);


%2D contour plot
figure(1);
contourf(Xq, Yq, Uq, 20, 'LineColor', 'none');
colorbar;
xlabel('x');
ylabel('y');
title('2D Poisson Solution Contour');
axis image;
