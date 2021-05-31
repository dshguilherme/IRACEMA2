clearvars
clc
close all
%% Example 3: Scordelis-Lo Roof

E = 4.32e8;
vu = 0.33;
g = @(x) [0 0 -90];
R = 25;
L = 50;
phi = 40*pi/180;
t = 0.25;
center = [0 0 0];
center2 = [0 L/2 0];
initial_point1 = [0 0 R-t/2];
initial_point2 = [0 0 R+t/2];
normal = [0 1 0];
arc1 = bs_arc(center,initial_point1,phi, normal);
arc2 = bs_arc(center,initial_point2,phi, normal);
section = bs_ruled_surface(arc1,arc2);
domain = bs_extrusion(section,[0 L/2 0]);

% Refinement
domain.degree_elevate(1,1);
domain.degree_elevate(1,2);
domain.degree_elevate(1,3);
Xi = linspace(0,1,16);
Xi = Xi(2:end-1);
domain.knot_refine(Xi,1);
domain.knot_refine(Xi,2);
domain.knot_refine(0.5,3);

% Assembly
asb = Elastic(E,vu,"gauss",3,domain);
K = asb.build_stiffness;
F = asb.force_vector(g);

% Boundary Conditions
id = asb.id_matrix;

P = domain.points;
P = cell2mat(P(:));

x_points = find(P(:,1) == 0);
x_dofs = id(x_points,1);

xz_points = find(P(:,2) == 0);
xz_dofs = id(xz_points,[1 3]);

clamped_dofs = unique([x_dofs(:); xz_dofs(:)]);

g = 0;
[d, F, solution] = asb.dirichlet_linear_solve(K,F,g,clamped_dofs);
solution.plot_solution(3)
colormap("hot");
