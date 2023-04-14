clearvars
close all
clc
% Geometry
P1 = [0 0 0 1];
P2 = [1 0 0 1];
U = [0 0 1 1];
line1 = Geometry(1, {U}, {P1, P2}, [1]);
line2 = bs_translation(line1,[0 1 0]);
domain = bs_ruled_surface(line1, line2);

% Material properties
lambda = 0.001;
mobility = 50;

% Refinement
Xi = linspace(0,1,5);
Xi = Xi(2:end-1);
domain.degree_elevate(1,1);
domain.degree_elevate(1,2);
domain.knot_refine(Xi,1);
domain.knot_refine(Xi,2);

% Assembler
max_steps = 1000;
theta = 1;
dt = 5e-6;
t_max = 5;
asb = CahnHilliardMixed(lambda, mobility, domain, theta, dt, t_max, max_steps);


%% Solving
asb.timeLoop;
asb.plotEvolution;
solution = TimeDependentSolution(asb, asb.solution_array);
solution.snapSolution('mixed_test',10);
  