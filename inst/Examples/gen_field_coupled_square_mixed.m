clearvars
close all
clc

%% Geometry
L = 1;
h = 1;
domain = bs_rectangle(L,h);

% Refinement
refinements = 4;
elevations = 1;
domain.uniform_k_refine(refinements, elevations);

%% Material properties
lambda = 0.001;
mobility = 1;
eta = 0.5;

% Steel
YOUNG = 200e9;
POISSON = 1/3;
rho = 8000;
alpha = 0;
beta = 0;

%% Assemblers
max_steps = 200;
theta = 1;
dt = 5e-6;
t_max = 5;
frequency = 0;
cf_asb = Elastodynamic(YOUNG, POISSON, rho, alpha, beta, "gauss", 2, domain);

% Boundary Conditions
b = domain.extract_boundaries;
clamped_cpoints = b{1,4};
id = cf_asb.id_matrix;
clamped_dofs = id(clamped_cpoints,:);
clamped_dofs = clamped_dofs(:);
g = zeros(size(clamped_dofs));
sides = [2];
trac = {[0 0 0]}; % Traction

h1 = 0.75;
h2 = 0.25;
intensity = 1e6;
P = domain.points;
dx = P{end} - P{end-1};
dx = dx(1);
force_L = L-dx;
force_function = @(x) cantileverForceAdjust(x, force_L, h1, h2, intensity);



gen_asb = GeneralizedPhasePhield(cf_asb, trac, ...
                sides, clamped_dofs, force_function, frequency, eta, lambda, ...
                mobility, domain, theta, dt, t_max, max_steps);

gen_asb.res_tol = 1e-3;
gen_asb.newton_max_steps = 40;
gen_asb.mode = "GMRES";
%% Solving
option = "elastic";
gen_asb.staggeredTimeLoop(option)
gen_asb.plotEvolution(option)
solution = TimeDependentSolution(gen_asb, gen_asb.solution_array);
clearvars -except gen_asb solution
save('gen_field_elastic.mat')
solution.snapSolution('gen_field_test2_',30);