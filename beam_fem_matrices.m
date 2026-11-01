function [K, M, freeDofs, x] = beam_fem_matrices(E, rho, A, I, L, nElem, boundaryCondition)
%BEAM_FEM_MATRICES Assemble Euler--Bernoulli beam FEM matrices.
% Nodal degrees of freedom are ordered as [w theta] at each node.

arguments
    E (1,1) double {mustBePositive}
    rho (1,1) double {mustBePositive}
    A (1,1) double {mustBePositive}
    I (1,1) double {mustBePositive}
    L (1,1) double {mustBePositive}
    nElem (1,1) double {mustBeInteger,mustBePositive}
    boundaryCondition (1,:) char = 'simply-supported'
end

nDof = 2 * (nElem + 1);
Le = L / nElem;
K = zeros(nDof);
M = zeros(nDof);
ke = E * I / Le^3 * [12 6*Le -12 6*Le; 6*Le 4*Le^2 -6*Le 2*Le^2; ...
    -12 -6*Le 12 -6*Le; 6*Le 2*Le^2 -6*Le 4*Le^2];
me = rho * A * Le / 420 * [156 22*Le 54 -13*Le; 22*Le 4*Le^2 13*Le -3*Le^2; ...
    54 13*Le 156 -22*Le; -13*Le -3*Le^2 -22*Le 4*Le^2];

for e = 1:nElem
    dofs = 2*e-1 : 2*e+2;
    K(dofs,dofs) = K(dofs,dofs) + ke;
    M(dofs,dofs) = M(dofs,dofs) + me;
end

switch lower(strtrim(boundaryCondition))
    case {'simply-supported', 'simply supported', 'ss'}
        constrained = [1, nDof-1];       % w(0) = w(L) = 0
    case {'cantilever', 'clamped-free', 'clamped free'}
        constrained = [1, 2];            % w(0) = theta(0) = 0
    otherwise
        error('Unsupported boundary condition: %s', boundaryCondition);
end
freeDofs = setdiff(1:nDof, constrained);
x = linspace(0, L, nElem+1).';
end
