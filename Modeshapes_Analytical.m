% Analytical modes and frequencies for the simply supported beam.
clearvars; close all; clc;
E = 210e9; rho = 7850; b = 0.05; h = 0.03; L = 1.0;
A = b*h; I = b*h^3/12; nModes = 6;
x = linspace(0,L,1000).';
omega = zeros(nModes,1);
figure('Name','Analytical mode shapes'); hold on;
for n = 1:nModes
    omega(n) = (n*pi/L)^2 * sqrt(E*I/(rho*A));
    plot(x, sin(n*pi*x/L), 'LineWidth', 1.1);
end
grid on; xlabel('x (m)'); ylabel('mode amplitude');
title('Simply supported analytical mode shapes');
legend(compose('Mode %d',1:nModes), 'Location','southwest');
disp(table((1:nModes).', omega/(2*pi), ...
    'VariableNames', {'Mode','Frequency_Hz'}));
