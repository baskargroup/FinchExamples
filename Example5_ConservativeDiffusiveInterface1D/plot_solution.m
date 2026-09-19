clc;close all; clear;
% Load the CSV file: x, phi (numerical), phi (exact)
data = readmatrix('solution.csv');
x       = data(:,1);
phi     = data(:,2);
phi_ex  = data(:,3);
% Plot both solutions
figure;
plot(x, phi, 'o-', 'DisplayName', 'Numerical'); hold on;
plot(x, phi_ex, 'r-', 'DisplayName', 'Exact (translated tanh)');
xlabel('x');
ylabel('\phi(x)');
title('1D CDI: interface advected by uniform velocity');
legend show;
grid on;
ylim([-0.1 1.1]);
fprintf('max error = %.4e\n', max(abs(phi - phi_ex)));
% Save the figure as PNG
saveas(gcf, 'cdi1d_solution.png')
