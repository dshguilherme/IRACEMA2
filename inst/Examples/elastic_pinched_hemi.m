clearvars
clc
close all
%% Geometry and Material
rho = 10.0;
t = 0.04;
E = 6.825e7;
vu = 0.3;

r = rho-t/2;
R = rho+t/2;

P{1,1,1} = [r 0 0 1]; P{1,1,2} = [R 0 0 1];
P{1,2,1} = [r 0 r 1/sqrt(2)]; P{1,2,2} = [R 0 R 1/sqrt(2)];
P{1,3,1} = [0 0 r 1]; P{1,3,2} = [0 0 R 1];
P{2,1,1} = [r r 0 1/sqrt(2)]; P{2,1,2} = [R R 0 1/sqrt(2)];
P{2,2,1} = [r r r 1/2]; P{2,2,2} = [R R R 1/2];
P{2,3,1} = [0 0 r 1/sqrt(2)]; P{2,3,2} = [0 0 R 1/sqrt(2)];
P{3,1,1} = [0 r 0 1]; P{3,1,2} = [0 R 0 1];
P{3,2,1} = [0 r r 1/sqrt(2)]; P{3,2,2} = [0 R R 1/sqrt(2)];
P{3,3,1} = [0 0 r 1]; P{3,3,2} = [0 0 R 1];

U = [0 0 0 1 1 1];
V = [0 0 0 1 1 1];
W = [0 0 1 1];
p = [2 2 1];

domain = Geometry(3,{U,V,W},P,p);

%% Refinement

% 2 quadratic elements in thickness:
domain.degree_elevate(1,3);
domain.knot_refine(0.5,3);  

% Quintic Elements on patches
domain.degree_elevate(3,1);
domain.degree_elevate(3,2);

% 15 elements per patch side
Xi = linspace(0,1,16);
Xi = Xi(2:end-1);
domain.knot_refine(Xi,1);
domain.knot_refine(Xi,2);

%% Assembly
asb = Elastic(E,vu,"gauss",3,domain);
K = asb.build_stiffness;
F = zeros(length(K),1);

%% Boundary Conditions
id = asb.id_matrix;
boundaries = domain.extract_boundaries;

% Symmetry on xz plane at y = 0
xz_points = boundaries{1,2};
xz_dofs = id(2,xz_points);

% Symmetry on yz plane at x = 0
yz_points = boundaries{2,2};
yz_dofs = id(1,yz_points);

% Fixed points on degenerate boundary
degen_points = boundaries{4,2};
degen_dofs = id(:,degen_points);

bcs = unique([xz_dofs(:); yz_dofs(:); degen_dofs(:)]);
g = zeros(numel(bcs),1);

% Point Forces at [R+t/2 0 0] and [0 R+t/2 0]
P = cell2mat(domain.points(:));
Fx = find(P(:,1) == R & P(:,2) == 0 & P(:,3) == 0);
Fx_dof = id(1,Fx);
Fy = find(P(:,1) == 0 & P(:,2) == R & P(:,3) == 0);
Fy_dof = id(2,Fy);

F(Fx_dof) = 0.5;
F(Fy_dof) = -0.5;
F = sparse(F);
[d, F, solution] = asb.dirichlet_linear_solve(K,F,g,bcs);
[x, d] = solution.eval_solution([0 0 1])
solution.plot_solution(1);
colormap("jet");

