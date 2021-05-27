clearvars
clc
%% Example 4: 2D Plate with circular hole

% Geometric and Material
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
Xi = linspace(0,1,36);
Xi = Xi(2:end-1);
domain.uniform_k_refine(Xi,2);
% domain.degree_elevate(4,1)
% domain.degree_elevate(4,2)
% domain.knot_refine([0.1 0.2 0.3 0.4 0.6 0.7 0.8 0.9],1)
% domain.knot_refine([0.1 0.2 0.3 0.4 0.6 0.7 0.8 0.9],2)

E = 1e5;
vu = 0.3;

% Assembly

asb = Elastic(E,vu,"gauss",2,domain);
K = asb.build_stiffness;

% Boundary Conditions
boundaries = domain.extract_boundaries;

x_sym_bc = boundaries(1,:);
y_sym_bc = boundaries(2,:);
stress_bc = boundaries(3,:);
zero_stress_bc = boundaries(4,:);

% Neumann BCs
id = asb.id_matrix;
[s1 s1] = size(K);
F = zeros(s1,1);
% h = [-1 0];
F = asb.constant_neumann_bc(F,h,stress_bc);
h = @(x) if x > -4 [0 0] else [-1 0] end;
F = asb.variable_neumann_bc(F,h,stress_bc);
% Dirichlet BCs

x_points = x_sym_bc{2};
x_dofs = id(x_points,1);

y_points = y_sym_bc{2};
y_dofs = id(y_points,2);

clamped_dofs = unique([x_dofs(:); y_dofs(:)]);
% clamped_dofs = x_dofs(:);
g = zeros(numel(clamped_dofs,1));

% Penalty to enforce the overlapping control points
% to have the same displacement
P = domain.points;
P = cell2mat(P(:));
[s1, ~] = size(domain.points);

ppoint = find(P(:,1) == -L & P(:,2) == L);
for i=ppoint(1):s1:length(P)
    pdofs = id([i i+1],:);
    penaltyStiffness = 1000*[1 -1; -1 1];
    px = pdofs(1,:);
    py = pdofs(2,:);
    K(px,px) = K(px,px) +penaltyStiffness;
    K(py,py) = K(py,py) +penaltyStiffness;
end

[d, F, solution] = asb.dirichlet_linear_solve(K,F,g,clamped_dofs);
solution.plot_solution(1)
[x, xx] = solution.eval_solution([0 1])
colormap(parula)
