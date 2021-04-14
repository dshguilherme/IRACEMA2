%% Example 6: Clamped rectangular plate with Elastodynamic PDE

% Geometric and material:

L = 1.5;
h = 1;
P1 = [0 0 0 1];
P2 = [L 0 0 1];
P3 = [0 h 0 1];
P4 = [L h 0 1];

U = [0 0 1 1];
V = U;
p = [1 1];

domain = Geometry(2, {U,V}, {P1,P2; {P3,P4},p);

domain = geo_extrusion(domain,[0 0 0.001]);

rho = 1;
E = 10e3;
vu = 0.3;

% Assembly

asb = Elastodynamic(rho, E, vu, "gauss", 3, domain);
[K, M] = asb.build_matrices;

% Boundary Conditions
[K_e, M_e] = asb.clamp_boundary(K,M);

% Solution
[d, omega] = eigs(K_e,M_e,'sm');