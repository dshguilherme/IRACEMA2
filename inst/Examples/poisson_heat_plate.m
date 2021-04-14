%% Example 1: Steady-state heat equation on a plate
% Neumann boundary condition at x = 0 (derivative = 0)
% Dirichlet boundary condition of 100C at y = 0
% The other two boundaries have a Robin Boundary condition r = hT/k, B = h/k
% h = 750W/m2C T = 0C (273.15K) k = 52W/mC
% The benchmark is temperature at [0.6,0.2] is 18.2C

% Geometric and material:
L = 0.6;
h = 1.0;
P1 = [0 0 0 1];
P2 = [L 0 0 1];
P3 = [0 h 0 1];
P4 = [L h 0 1];

U = [0 0 1 1];
V = [0 0 1 1];
p = [1 1];

domain = Geometry(2, {U,V}, {P1,P2; P3,P4},p);

alpha = 52;

% Refinement
% domain.k_refine;

% Assembly
asb = Poisson(alpha,"gauss",1,domain);
K = asb.build_stiffness;
F = zeros(length(K),1);
F = sparse(F);
% Boundary conditions
boundaries = domain.extract_boundaries;

b_robin1 = boundaries(3,:);
b_robin2 = boundaries(4,:);
b_robin = [b_robin1; b_robin2];
r = 0;
beta = 750/alpha;

[K,F] = asb.robin_bc(K,F,r,beta,b_robin);

% Dirichlet and Solution
b_dirichlet = cell2mat(boundaries(2,2));
g = 100;

[d, ~, solution] = asb.dirichlet_linear_solve(K,F,g,b_dirichlet);

%Post processing
solution.plot_solution;
point = solution.eval_solution(1,0.2)
