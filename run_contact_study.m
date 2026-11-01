%RUN_CONTACT_STUDY Sweep contact stiffness and plot hysteresis/energy.
% Part 1 remains simply supported. Part 2 uses a fixed-free cantilever mode.

clearvars; close all; clc;
E = 210e9; rho = 7850; b = 0.05; h = 0.03; L = 1;
 A = b*h; I = b*h^3/12;
beam = cantilever_modal_properties(E,rho,A,I,L);
 mModal = beam.modalMass;
omega1 = beam.naturalFrequency;

base = struct('mass',mModal,'damping',0.02*mModal*omega1, ...
    'naturalFrequency',omega1,'forceAmplitude',5.0,'driveFrequency',0.9*omega1, ...
    'duration',12/0.9/omega1*2*pi,'timeStep',1e-4,'gap',1e-6, ...
    'frictionCoefficient',0.35,'normalLoad',mModal*9.81, ...
    'contactLocation',L,'modalShapeAtContact',beam.modeShapeAtTip, ...
    'frictionVelocity',1e-5,'youngsModulus',E,'beamLength',L, ...
    'beamThickness',h,'curvatureAtStressLocation', ...
    beam.curvatureAtRootPerDisplacement);
stiffnesses = [1e3 1e4 1e5 1e6];
results = cell(size(stiffnesses));
for k = 1:numel(stiffnesses)
    base.contactStiffness = stiffnesses(k);
    results{k} = frictional_contact_response(base);
end

figure('Name','Contact hysteresis'); hold on;
for k = 1:numel(results)
    r = results{k};
    plot(r.contactDisplacement, r.contactForce, ...
        'DisplayName',sprintf('k_c = %.0e N/m',stiffnesses(k)));
end
grid on; xlabel('contact displacement (m)'); ylabel('contact force (N)');
title('Frictional contact force-displacement hysteresis');
legend('Location','best');

energy = cellfun(@(r) r.energyDissipated,results);
fprintf('Contact stiffness (N/m)    dissipated energy (J)\n');
disp([stiffnesses(:), energy(:)]);
figure('Name','Contact stiffness sweep');
semilogx(stiffnesses,energy,'-o','LineWidth',1.2); grid on;
xlabel('contact stiffness (N/m)'); ylabel('dissipated energy (J)');
title('Energy dissipation versus contact stiffness');
