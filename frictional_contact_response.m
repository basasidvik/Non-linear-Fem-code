function result = frictional_contact_response(p)
%FRICTIONAL_CONTACT_RESPONSE Simulate a beam mode with frictional contact.
% Contact is unilateral: it activates only when the contact-point
% displacement exceeds p.gap. The normal force is elastic (p.contactStiffness)
% and the tangential friction force opposes velocity, producing a closed
% force-displacement hysteresis loop.

required = {'mass','damping','naturalFrequency','forceAmplitude', ...
    'driveFrequency','duration','timeStep','contactStiffness'};
for i = 1:numel(required)
    if ~isfield(p,required{i})
        error('Missing parameter: %s', required{i});
    end
end
defaults = struct('gap',0,'frictionCoefficient',0.35,'normalLoad',0, ...
    'contactLocation',0.5,'modalShapeAtContact',1,'frictionVelocity',1e-5, ...
    'youngsModulus',210e9,'beamLength',1,'beamThickness',0.03, ...
    'curvatureAtStressLocation',-(pi)^2);
names = fieldnames(defaults);
for i = 1:numel(names)
    if ~isfield(p,names{i}), p.(names{i}) = defaults.(names{i}); end
end
validateattributes(p.timeStep, {'numeric'}, {'scalar','positive'});
validateattributes(p.duration, {'numeric'}, {'scalar','positive'});

if p.normalLoad == 0
    p.normalLoad = p.mass * 9.81;
end
time = (0:p.timeStep:p.duration).';
shape = p.modalShapeAtContact;
state0 = [0; 0];
rhs = @(t,y) contact_rhs(t,y,p,shape);
options = odeset('RelTol',1e-7,'AbsTol',1e-10);
[~, state] = ode45(rhs,time,state0,options);
q = state(:,1); qdot = state(:,2);
contactDisplacement = shape*q;
contactVelocity = shape*qdot;
penetration = max(contactDisplacement-p.gap,0);
normalForce = p.contactStiffness*penetration;
frictionForce = zeros(size(time));
inContact = penetration > 0;
frictionForce(inContact) = -p.frictionCoefficient*p.normalLoad .* ...
    tanh(contactVelocity(inContact)/p.frictionVelocity);
contactForce = normalForce + frictionForce;

result.time = time;
result.modalDisplacement = q;
result.modalVelocity = qdot;
result.contactDisplacement = contactDisplacement;
result.contactForce = contactForce;
result.normalForce = normalForce;
result.frictionForce = frictionForce;
result.penetration = penetration;
result.inContact = inContact;
result.energyDissipated = trapz(time, max(-frictionForce.*contactVelocity,0));
result.hysteresisWork = abs(trapz(contactDisplacement,contactForce));
result.contactStiffness = p.contactStiffness;
% Bending stress/strain at the requested stress location. The caller can
% provide the cantilever root curvature, which is more useful than the zero
% curvature at a mathematically free tip.
curvature = p.curvatureAtStressLocation*q;
result.bendingStrain = -0.5*p.beamThickness*curvature;
result.bendingStress = p.youngsModulus*result.bendingStrain;
end

function dydt = contact_rhs(t,y,p,shape)
q = y(1); qdot = y(2);
x = shape*q; v = shape*qdot;
penetration = max(x-p.gap,0);
normal = p.contactStiffness*penetration;
friction = 0;
if penetration > 0
    friction = -p.frictionCoefficient*p.normalLoad*tanh(v/p.frictionVelocity);
end
contactGeneralized = shape*(normal+friction);
forcing = p.forceAmplitude*sin(p.driveFrequency*t)*shape;
dydt = [qdot; (forcing-p.damping*qdot-p.mass*p.naturalFrequency^2*q- ...
    contactGeneralized)/p.mass];
end
