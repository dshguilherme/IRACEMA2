clearvars
clc
%% Example 4: 2D Plate with circular hole

%    stress
% s ____________
% t|            |
% r|            | dx = 0
% e|          __|
% s|_________|
% s   dy = 0
% 

% Geometry and Material
L = 4;
P1 = [-L 0 0 1];
P2 = [-L L 0 1];
P3 = [0 L 0 1];
U = [0 0 0.5 1 1];

line = Geometry(1, {U}, {P3,P2,P1}, [1]);
R = 1;
center = [0 0 0];
initial_point = [0 R 0];
theta = pi/2;
normal = [0 0 1];
arc = bs_arc(center, initial_point, theta, normal);

domain = bs_ruled_surface(line,arc);

% Refinement
refinements = 10;
elevations = 0;
domain.uniform_k_refine(refinements, elevations);

% Material Properties
YOUNG = 1e5;
POISSON = 0.3;
rho = 1;
alpha = 0;
beta = 0;

%% Assembler
asb = Elastodynamic(YOUNG, POISSON, rho, alpha, beta, "gauss", 2, domain);

%% Boundary Conditions
b = domain.extract_boundaries;

x_clamp_points = b{1,4};
y_clamp_points = b{2,4};

id = asb.id_matrix;

x_clamp_dofs = id(x_clamp_points,1);
y_clamp_dofs = id(y_clamp_points,2);

clamped_dofs = [x_clamp_dofs(:); y_clamp_dofs(:)];

% Traction
Tx = 10;
sides = [3];
t = {@(x) exactStress(Tx,R,x)};

%% Solving
d = asb.solveLinearElasticity(@(x) [0 0 0], t, sides, 0, clamped_dofs);

% Solution
sol = Solution(asb,d);
L2 = asb.errorL2(d, @(x) exactStress(Tx,R,x),2)
