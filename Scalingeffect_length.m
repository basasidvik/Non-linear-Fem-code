% Length scaling of simply supported analytical frequencies.
clearvars; close all; clc;
E = 210e9; rho = 7850; b = 0.05; h = 0.03;
A = b*h; I = b*h^3/12; lengths = 1:0.5:2; nModes = 6;
figure('Name','Length scaling'); hold on;
for L = lengths
    n = (1:nModes).';
    omega = (n*pi/L).^2 * sqrt(E*I/(rho*A));
    plot(n, omega, '-o', 'DisplayName', sprintf('L = %.1f m', L));
end
grid on; xlabel('mode number'); ylabel('natural frequency (rad/s)');
title('Simply supported length scaling'); legend('Location','northwest');
