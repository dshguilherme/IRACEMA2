clearvars
close all
clc
%% Geometry
L = pi;
h = exp(1);
P1 = [0 0 0 1];
P2 = [L 0 0 1];
U = [0 0 1 1];
line1 = Geometry(1, {U}, {P1, P2}, [1]);
line2 = bs_translation(line1,[0 h 0]);
domain = bs_ruled_surface(line1, line2);

% Refinement
Xi = linspace(0,1,30);
Xi = Xi(2:end-1);
domain.degree_elevate(1,1);
domain.degree_elevate(1,2);
domain.knot_refine(Xi,1);
domain.knot_refine(Xi,2);

% Material properties for Steel
YOUNG = 200e9;
POISSON = 1/3;
rho = 8000;
alpha = 0;
beta = 0;



%% Assembler
asb = Elastodynamic(YOUNG, POISSON, rho, alpha, beta, "gauss", 2, domain);

%% Boundary Conditions
b = domain.extract_boundaries;
clamped_cpoints = b{1,4};
id = asb.id_matrix;
clamped_dofs = id(clamped_cpoints,:);
clamped_dofs = clamped_dofs(:);
g = zeros(size(clamped_dofs));

t = {[0 0 0]}; % Traction
sides = [2];

%% Solutions

% Elasticity
d_traction = asb.solveLinearElasticity(@(x) [0 0 0], t, sides, g, clamped_dofs);
d_force = asb.solveLinearElasticity(@(x) cantileverForce(x), t, sides, g, clamped_dofs);
[K, M, F] = asb.assembleSystem(@(x) cantileverForce(x));
% Modal
num_modes = 10;
[modes, omega, modal_solution] = asb.naturalModes(sides, num_modes);

% Dynamic
frequency = 6000; 
d_dynamic = asb.solveDynamicProblem(@(x) cantileverForce(x), frequency, clamped_dofs);

frequency_range = 500:100:10000;
points = b{2,4};
FRFs = asb.pointFRF(@(x) cantileverForce(x), frequency_range, clamped_dofs, points);


