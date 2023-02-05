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
lambda = 0.05;
theta = 3/2;
diffusivity = 50;

% Refinement
Xi = linspace(0,1,10);
Xi = Xi(2:end-1);
domain.degree_elevate(1,1);
domain.degree_elevate(1,2);
domain.knot_refine(Xi,1);
domain.knot_refine(Xi,2);

% Initial condition
initial_c = 0.5 +0.02.*randn(numel(domain.points),1);

% Assembler
asb = CahnHilliard(lambda,theta,diffusivity,initial_c,domain);

% Solving
c = asb.solveSystem(2e-6,150);
c = c(:,any(c));
timed_solutions = TimeDependentSolution(asb,c);
F = timed_solutions.recordSolution("ch_square.mp4",2);