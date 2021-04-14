%% Example 2: Modes of vibration of a beam
% Dirichlet boundary conditions at x = 0 and x = L (clamped)

% Geometric and material:

L = 1;
P1 = [0 0 0 1];
P2 = [L 0 0 1];
U = [0 0 1 1];
domain = Geometry(1, {U}, {P1,P2}, 1);

alpha = 1;
density = 1;

% Refinement
% domain.k_refine

% Assembly
asb = Helmholtz(density,alpha,"gauss",1,domain);
[K,M] = asb.build_matrices;

% Boundary Conditions
boundaries = domain.extract_boundaries;
[K_c, M_c] = asb.clamp_boundary(K,boundaries);

% Solution
[d, omega] = eigs(K,M,'sm',100);
omega = sqrt(diag(omega));

%Post processing