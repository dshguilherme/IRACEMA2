clearvars
clc
close all
%% Example 3: Scordelis-Lo Roof

E = 4.32e8;
vu = 0.33;
g = [0 0 -90];
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
domain.degree_elevate(3,1);
domain.degree_elevate(3,2);
domain.degree_elevate(1,3);
Xi = linspace(0,1,16);
Xi = Xi(2:end-1);
domain.knot_refine(Xi,1);
domain.knot_refine(Xi,2);
domain.knot_refine(0.5,3);

% Assembly
asb = Elastic(E,vu,"gauss",3,domain);
K = asb.build_stiffness;
F = asb.constant_force(g);

% Boundary Conditions
id = asb.id_matrix;
boundaries = domain.extract_boundaries;

rigid_d = boundaries(5,:);
rigid_points = rigid_d{2};
rigid_dofs = id([1 3],rigid_points);

sym_plane1 = boundaries(1,:);
x_points = sym_plane1{2};
x_dofs = id(1,x_points);

sym_plane2 = boundaries(6,:);
y_points = sym_plane2{2};
y_dofs = id(2,y_points);

clamped_dofs = unique([rigid_dofs(:); x_dofs(:); y_dofs(:)]);
bc = zeros(numel(clamped_dofs),1);
[d, F, solution] = asb.dirichlet_linear_solve(K,F,bc,clamped_dofs);
solution.plot_solution(3)
colormap("hot");
