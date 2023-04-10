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
lambda = 0.0005;
mobility = 1;

% Refinement
Xi = linspace(0,1,10);
Xi = Xi(2:end-1);
domain.degree_elevate(1,1);
domain.degree_elevate(1,2);
domain.knot_refine(Xi,1);
domain.knot_refine(Xi,2);

% Initial condition
initial_c = 0.5 +0.05.*randn(numel(domain.points),1);

% Assembler
max_steps = 1000;
theta = 1;
dt = 5e-6;
t_max = 1;
asb = CahnHilliardMixed(lambda, mobility, domain, theta, dt, t_max, max_steps);


%% Solving
[~ , ~] = asb.timeLoop;
solution = TimeDependentSolution(asb, asb.solution_array);
F = solution.recordSolution('mixed_test',30);
  