function beam = cantilever_modal_properties(E, rho, A, I, L, nPoints)
%CANTILEVER_MODAL_PROPERTIES First bending mode of a fixed-free beam.
% The mode is normalized to unit displacement at the free tip.
if nargin < 6, nPoints = 2001; end
betaL = 1.875104068711961;
beta = betaL/L;
x = linspace(0,L,nPoints).';
sigma = (cosh(betaL)+cos(betaL))/(sinh(betaL)+sin(betaL));
phi = cosh(beta*x)-cos(beta*x)-sigma*(sinh(beta*x)-sin(beta*x));
phi = phi/phi(end);
modalMass = trapz(x,rho*A*phi.^2);
beam.x = x;
beam.modeShape = phi;
beam.modeShapeAtTip = 1;
beam.modalMass = modalMass;
beam.naturalFrequency = beta^2*sqrt(E*I/(rho*A));
beam.beta = beta;
% Curvature at the fixed root, used for bending stress/strain reporting.
beam.curvatureAtRootPerDisplacement = 2*beta^2;
end
