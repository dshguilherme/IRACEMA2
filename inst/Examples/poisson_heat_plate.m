%% Example 2: Steady-state heat equation on a plate
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
U = [0 0 1 1];
line1 = Geometry(1, {U}, {P1,P2}, [1]);
line2 = bs_translation(line1,[0 h 0]);
domain = bs_ruled_surface(line1,line2);
alpha = 52;
h = 750;
T_inf = 273.15;

% Refinement
Xi = linspace(0,1,16);
Xi = Xi(2:end-1);
domain.uniform_k_refine(Xi,1);

% Assembly
asb = Poisson(alpha,"gauss",1,domain);
K = asb.build_stiffness;
F = zeros(length(K),1);
F = sparse(F);
% Boundary conditions
boundaries = domain.extract_boundaries;

b_robin1 = boundaries(2,:);
b_robin2 = boundaries(4,:);
b_robin = [b_robin1; b_robin2];
r = h*T_inf;
beta = h;

[K,F] = asb.robin_bc(K,F,r,beta,b_robin);

% Dirichlet and Solution
b_dirichlet = cell2mat(boundaries(3,2));
g = 373.15;

[d, ~, solution] = asb.dirichlet_linear_solve(K,F,g,b_dirichlet);

%Post processing
solution.plot_solution(1);
colormap("hot");
xlim([-.2 0.8]);
pbaspect([1 1 1]);
[x, d] = solution.eval_solution([1,0.2]);
d-273.15
