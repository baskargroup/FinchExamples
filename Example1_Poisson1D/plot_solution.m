clc;close all; clear;

% Load the CSV file
data = readmatrix('solution.csv');

x = data(:,1);         % x coordinates
u = data(:,2);         % numerical solution

% Analytical solution
u_exact = sin(10*pi*x) .* sin(pi*x);

% Plot both
% Plot both solutions
figure;
plot(x, u, 'o-', 'DisplayName', 'Numerical'); hold on;
plot(x, u_exact, 'r-', 'DisplayName', 'Analytical');
xlabel('x');
ylabel('u(x)');
title('1D Poisson Solution vs Analytical');
legend show;
grid on;

% Save the figure as PNG
saveas(gcf, 'poisson1d_solution.png')  

% Or save as PDF
% saveas(gcf, 'poisson1d_solution.pdf')  
