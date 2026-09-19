clc;clear;close all
% Final state of the 2D CDI solution.
% solution.csv: x, y, phi (numerical), phi (exact)
data  = readmatrix('solution.csv');
x     = data(:,1);
y     = data(:,2);
phi   = data(:,3);
phi_ex = data(:,4);

% Interpolate onto a regular grid for plotting
nx = 200; ny = 200;
xq = linspace(min(x), max(x), nx);
yq = linspace(min(y), max(y), ny);
[Xq, Yq] = meshgrid(xq, yq);
F  = scatteredInterpolant(x, y, phi, 'linear', 'none');
Fe = scatteredInterpolant(x, y, phi_ex, 'linear', 'none');
PHI  = F(Xq, Yq);
PHIe = Fe(Xq, Yq);

% Numerical field with phi = 0.5 contours (numerical vs exact)
figure(1);
contourf(Xq, Yq, PHI, 20, 'LineColor', 'none'); hold on;
contour(Xq, Yq, PHI,  [0.5 0.5], 'k-', 'LineWidth', 1.5);
contour(Xq, Yq, PHIe, [0.5 0.5], 'r--', 'LineWidth', 1.5);
colorbar; caxis([0 1]);
xlabel('x'); ylabel('y');
title('2D CDI: final \phi (black: numerical, red dashed: exact interface)');
axis image;
saveas(gcf, 'cdi2d_solution.png');

% Pointwise error
figure(2);
contourf(Xq, Yq, abs(PHI - PHIe), 20, 'LineColor', 'none');
colorbar;
xlabel('x'); ylabel('y');
title('Absolute error');
axis image;

fprintf('max nodal error = %.4e\n', max(abs(phi - phi_ex)));
