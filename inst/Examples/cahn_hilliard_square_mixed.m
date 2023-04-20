clearvars
close all
clc

% Geometry
L = 1;
h = 1;
domain = bs_rectangle(L,h);

% Refinement
refinements = 8;
elevations = 1;
domain.uniform_k_refine(refinements, elevations);

% Material properties
lambda = 0.001;
mobility = 1;

% Assembler
max_steps = 3000;
theta = 1;
dt = 5e-6;
t_max = 5;
asb = CahnHilliardMixed(lambda, mobility, domain, theta, dt, t_max, max_steps);
asb.res_tol = 1e-3;
asb.newton_max_steps = 40;

%% Solving
asb.timeLoop;
asb.plotEvolution;
solution = TimeDependentSolution(asb, asb.solution_array);
clearvars -except asb solution
save('mixed_test.mat')
solution.snapSolution('mixed_test',10);
  