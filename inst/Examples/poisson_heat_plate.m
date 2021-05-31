%% Example 2: Steady-state heat equation on a plate
% Neumann boundary condition at x = 0 (derivative = 0)
% Dirichlet boundary condition of 100C at y = 0
% The other two boundaries have a Robin Boundary condition r = hT/k, B = h/k
% h = 750W/m2C T = 0C (273.15K) k = 52W/mC
% The benchmark is temperature at [0.6,0.2] is 18.2C
clearvars
clc
close all
% Geometric and material:
L = 0.6;
h = 1.0;
P1 = [0 0 0 1];
P2 = [L 0 0 1];
U = [0 0 1 1];
line1 = Geometry(1, {U}, {P1,P2},[1]);
line2 = bs_translation(line1,[0 h 0]);
domain = bs_ruled_surface(line1,line2);
alpha = 52;
h = 750;
T_inf = 0;
T_d = 100;

% Refinement
domain.degree_elevate(1,1);
domain.degree_elevate(1,2);
Xi = linspace(0,1,15);
Xi = Xi(2:end-1);
domain.knot_refine(Xi,1);
Xi = linspace(0,1,9);
Xi = Xi(2:end-1);
domain.knot_refine(Xi,2);
% Assembly
asb = Poisson(alpha,"gauss",1,domain);
K = asb.build_stiffness;
% Boundary conditions
boundaries = domain.extract_boundaries;

b_robin1 = boundaries(2,:);
b_robin2 = boundaries(4,:);
b_robin = [b_robin1; b_robin2];
r = h*T_inf;
beta = h;

[Kr,F] = asb.robin_bc(r,beta,b_robin);
K = K+Kr;

% Dirichlet and Solution
id = asb.id_matrix;
P = domain.points;
P = cell2mat(P(:));
b_points = find(P(:,2) == 0);
b_dirichlet = id(b_points);
g = T_d;

[d, ~, solution] = asb.dirichlet_linear_solve(K,F,g,b_dirichlet);

%Post processing
solution.plot_solution(1);
colormap("jet");
xlim([-.2 0.8]);
ylim([-.1 1.1]);
pbaspect([1 1 1]);
[x, d] = solution.eval_solution([1,0.4])

