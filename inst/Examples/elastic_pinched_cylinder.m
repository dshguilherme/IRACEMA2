clearvars
clc
clear all

%% Geometry and Material

R = 300;
P = 1;
L = 600;
E = 3e6;
vu = 0.3;
t = 3;

center = [0 0 0];
initial_point = [R+t/2 0 0];
line1 = bs_arc(center, initial_point, -pi/2, [0 1 0]);
line2 = bs_arc(center, [R-t/2 0 0], -pi/2, [0 1 0]);
cross_sec = bs_ruled_surface(line2, line1);

domain = bs_extrusion(cross_sec, [0 L/2 0]);

%% Refinement
domain.degree_elevate(3,1);
domain.degree_elevate(1,2);
domain.degree_elevate(4,3);

domain.knot_refine([0.5],2);
Xi = linspace(0,1,21);
Xi = Xi(2:end-1);
domain.knot_refine(Xi,1);
domain.knot_refine(Xi,3);

%% Assembly
asb = Elastic(E,vu,"gauss",3,domain);

K = asb.build_stiffness;
F = zeros(length(K),1);

%% Boundary Conditions
id = asb.id_matrix;

boundaries = domain.extract_boundaries;

% Rigid diaphragm at y = 0
diaphragm = boundaries(5,:);
dph_points = cell2mat(diaphragm(2));
dph_dofs = id([1 3],dph_points); % xz plane constrained

% Symmetry at z = 0
xy_sym = boundaries(1,:);
xy_points = cell2mat(xy_sym(2));
xy_dofs = id([3],xy_points);

% Symmetry at x = 0
yz_sym = boundaries(2,:);
yz_points = cell2mat(yz_sym(2));
yz_dofs = id([1],yz_points);

bcs = unique([dph_dofs(:); xy_dofs(:); yz_dofs(:)]);
g = zeros(numel(bcs),1);
% Point Force at [0 L R+t/2]
F_point = numel(domain.points); % It's the last cpoint
F_dof = id(3,F_point);
F(F_dof) = -P/4;
F = sparse(F);

%% Solution
[d, F, solution] = asb.dirichlet_linear_solve(K,F,g,bcs);
[x, d] = solution.eval_solution([1 1 1])



