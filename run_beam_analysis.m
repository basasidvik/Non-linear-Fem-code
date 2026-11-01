%RUN_BEAM_ANALYSIS FEM validation and nonlinear modal dynamics.
% The default case is the requested simply supported Euler--Bernoulli beam.
% The nonlinear transient is a reduced-order model based on mid-plane
% stretching (von Karman axial strain).

clearvars; close all; clc;

E = 210e9; rho = 7850; b = 0.05; h = 0.03; L = 1.0;
A = b*h; I = b*h^3/12; nElem = 40; nModes = 6;
boundaryCondition = 'simply-supported';

[K, M, freeDofs, x] = beam_fem_matrices(E, rho, A, I, L, nElem, boundaryCondition);
[V, D] = eig(K(freeDofs,freeDofs), M(freeDofs,freeDofs));
[lambda, order] = sort(real(diag(D)), 'ascend');
V = V(:,order);
omegaFEM = sqrt(lambda(1:nModes));

Phi = zeros(size(K,1), nModes);
for r = 1:nModes
    mode = V(:,r);
    mode = mode / sqrt(mode' * M(freeDofs,freeDofs) * mode);
    fullMode = zeros(size(K,1),1);
    fullMode(freeDofs) = mode;
    [~, peak] = max(abs(fullMode(1:2:end)));
    fullMode = sign(fullMode(2*peak-1)) * fullMode;
    Phi(:,r) = fullMode;
end

modeNumber = (1:nModes).';
omegaExact = (modeNumber*pi/L).^2 * sqrt(E*I/(rho*A));
fprintf('Simply supported beam: %d elements, %d modes\n', nElem, nModes);
fprintf('Mode       FEM (Hz)       analytical (Hz)      error (%%)\n');
fprintf('%4d %15.6f %20.6f %13.3f\n', [(1:nModes).', omegaFEM/(2*pi), ...
    omegaExact/(2*pi), 100*(omegaFEM-omegaExact)./omegaExact].');

figure('Name','FEM mode shapes');
plot(x, Phi(1:2:end,1:nModes), 'LineWidth', 1.1); grid on;
xlabel('x (m)'); ylabel('unit-mass mode amplitude');
title('Simply supported beam: FEM mode shapes');
legend(compose('Mode %d', 1:nModes), 'Location','southwest');

forceAmplitude = 1.0;                  % N, harmonic point force
driveFrequency = 0.80 * omegaExact(1); % rad/s
tEnd = 2.0; dt = 1e-3;
[time, q] = nonlinear_modal_response(E, rho, A, L, omegaExact, nModes, ...
    forceAmplitude, driveFrequency, L/2, tEnd, dt);

figure('Name','Nonlinear modal response');
plot(time, q, 'LineWidth', 1.0); grid on;
xlabel('time (s)'); ylabel('modal displacement (m)');
title('Geometrically nonlinear reduced-order response');
legend(compose('q_%d', 1:nModes), 'Location','best');

% Frictional contact is a separate, reusable extension of this beam model.
% Run run_contact_study.m to sweep contact stiffness and obtain hysteresis.

function [time, q] = nonlinear_modal_response(E, rho, A, L, omega, nModes, F0, ...
        driveFrequency, loadLocation, tEnd, dt)
%NONLINEAR_MODAL_RESPONSE Integrate multimode stretching nonlinearities.

xq = linspace(0, L, 4001).';
dphi = zeros(numel(xq), nModes);
for i = 1:nModes
    dphi(:,i) = (i*pi/L) * cos(i*pi*xq/L);
end
B = zeros(nModes);
for i = 1:nModes
    for j = 1:nModes
        B(i,j) = trapz(xq, dphi(:,i).*dphi(:,j));
    end
end
modalMass = rho * A * L / 2;
loadShape = sin((1:nModes).' * pi * loadLocation/L);
initial = zeros(2*nModes,1);
rhs = @(t,y) [y(nModes+1:end); -omega.^2 .* y(1:nModes) ...
    - (E*A/(2*L*modalMass)) * (y(1:nModes)'*B*y(1:nModes)) * B*y(1:nModes) ...
    + (F0/modalMass)*loadShape*cos(driveFrequency*t)];
time = (0:dt:tEnd).';
options = odeset('RelTol',1e-7, 'AbsTol',1e-9);
[~, y] = ode45(rhs, time, initial, options);
q = y(:,1:nModes);
end
