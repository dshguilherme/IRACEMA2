%% Example 5: Thin Beam with Elastodynamic PDE

% Geometric and material:

L = 1;
h = 1;
P1 = [0 0 0 1];
P2 = [L 0 0 1];
P3 = [0 h 0 1];
P4 = [L h 0 1];

U = [0 0 1 1];
V = U;
p = [1 1];

domain = Geometry(2, {U,V}, {P1,P2; {P3,P4},p);

domain = geo_extrusion(domain,[0 0 20*L]);

rho = 1;
E = 10e3;
vu = 0.3;

% Assembly

asb = Elastodynamic(rho, E, vu, "gauss", 3, domain);
[K, M] = asb.build_matrices;

% Boundary Conditions
boundaries = domain.extract_boundaries;
dirichlet_boundary = boundaries(1, 2);
dirichlet_boundary = unique(dirichlet_boundary);
[K_e, M_e] = asb.clamp_dofs(K,M,dirichlet_boundary);

% Solution
[d, omega] = eigs(K_e,M_e,'sm');