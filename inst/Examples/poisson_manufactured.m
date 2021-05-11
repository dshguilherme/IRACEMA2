%% Example 1: Poisson with method of manufactured solutions

% Solution: u = 1 + x^2  +2y^2
% -lap(u) = f -> f = -6
% with Dirichlet Conditions = 1 + x^2 +2y^2

% Geometry
P1 = [0 0 0 1];
P2 = [1 0 0 1];
U = [0 0 1 1];
line1 = Geometry(1, {U}, {P1, P2}, [1]);
line2 = bs_translation(line1,[0 1 0]);

domain = bs_ruled_surface(line1, line2);

% Refinement
Xi = linspace(0,1,11);
Xi = Xi(2:end-1);
domain.uniform_k_refine(Xi,1);

% Assembly
asb = Poisson(1,"gauss",1,domain);
K = asb.build_stiffness;
f = -6;
F = asb.constant_force(f);

% Boundary Conditions
boundaries = domain.extract_boundaries;
cpoints = boundaries(:,2);
cpoints = cell2mat(cpoints(:));
cpoints = unique(cpoints);
P = domain.points;
dirichlet_cp = P(cpoints);
u = @(x) 1 +(x(1).^2) +2*(x(2).^2);
g = cellfun(u,dirichlet_cp);
[d, F, solution] = asb.dirichlet_linear_solve(K,F,g,cpoints);
