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
Xi = linspace(0,1,18);
Xi = Xi(2:end-1);
domain.uniform_k_refine(Xi,2);
% domain.knot_refine(Xi,1)
% domain.knot_refine(Xi,2)

E = 1e5;
vu = 0.3;

% Assembly

asb = Elastic(E,vu,"gauss",2,domain);
K = asb.build_stiffness;

% Boundary Conditions
boundaries = domain.extract_boundaries;
stress_bc = boundaries(3,:);

% Neumann BCs
id = asb.id_matrix;
[s1 s1] = size(K);
F = zeros(s1,1);
Tx = 10;
h = @(x) plate_stress(x,Tx);
% h = @(x) [-Tx 0];
F = asb.neumann_bc(h,stress_bc);
% Dirichlet BCs

P = domain.points;
P = cell2mat(P(:));

x_points = find(P(:,1) == 0);
x_dofs = id(x_points,1);

y_points = find(P(:,2) == 0);
y_dofs = id(y_points,2);

clamped_dofs = [x_dofs(:); y_dofs(:)];

ppoint = find(P(:,1) == -L & P(:,2) == L);
for i=ppoint(1)%:s1:length(P)
    pdofs = id([i i+1],:);
    penaltyStiffness = 1000*[1 -1; -1 1];
    px = pdofs(1,:);
    py = pdofs(2,:);
    K(px,px) = K(px,px) +penaltyStiffness;
    K(py,py) = K(py,py) +penaltyStiffness;
end

g = 0;
[d, F, solution] = asb.dirichlet_linear_solve(K,F,g,clamped_dofs);
figure(1)
solution.plot_solution(1)
[x, xx] = solution.eval_solution([0 1]);
colormap(jet)
% figure(2)
% [stress_d, s_sol] = asb.project_stress(d);
% s_sol.plot_solution(1);



