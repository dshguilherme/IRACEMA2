clearvars
clc
%% Example 7: 3D Beam via Elastodynamic Equations
% Material Constants
rho = 100; % Steel beam
E = 1e3;
vu = 0.33;

% Geometry
L = 20;
b = 1;

P1 = [0 0 0];
P2 = [b 0 0];

line = bs_line(P1,P2);
line2 = bs_translation(line,[0 0 b]);

section = bs_ruled_surface(line,line2);
domain = bs_extrusion(section,[0 L 0]);

% Refinement
domain.degree_elevate(1,1);
domain.degree_elevate(1,2);
domain.degree_elevate(1,3);

domain.knot_refine([1/3 2/3],1);
domain.knot_refine([1/3 2/3],2);

Xi = linspace(0,1,101);
Xi = Xi(2:end-1);
domain.knot_refine(Xi,3);

% Assembly
asb = Elastodynamic(rho,E,vu,"gauss",3,domain);
[K, M] = asb.build_matrices;

% % BCs
id = asb.id_matrix;
boundaries = domain.extract_boundaries;
points = boundaries{5,2};
clamp_dofs = id(:,points);

% Solution
n_modes = 10;
[~, omega, solution] = asb.eigensolve(K,M,n_modes,clamp_dofs);

