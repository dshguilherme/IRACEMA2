%% Example 3: Eigenvalues of a membrane
% Clamped at the edges

% Geometric and material:
L = 1.5;
h = 1.0;
P1 = [0 0 0 1];
P2 = [L 0 0 1];
P3 = [0 h 0  1];
P4 = [L h 0 1];

U = [0 0 1 1];
V = [0 0 1 1];
p = [1 1];

domain = Geometry(1, {U,V}, {P1,P2; P3,P4},p);

alpha = 1;
density = 1;

% Refinement
domain.k_refine

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


