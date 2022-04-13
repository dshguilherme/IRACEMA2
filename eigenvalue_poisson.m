clearvars
close all
clc

% Geometry 
    P1 = [0 0 0 1];
    P2 = [1.5 0 0 1];
    U = [0 0 1 1];
    line1 = Geometry(1, {U}, {P1,P2}, [1]);
    line2 = bs_translation(line1, [0 1 0]);
    domain = bs_ruled_surface(line1, line2);

% Refinement
    Xi = linspace(0,1,10);
    Xi = Xi(2:end-1);
    domain.degree_elevate(1,1);
    domain.degree_elevate(1,2);
    domain.knot_refine(Xi,1);
    domain.knot_refine(Xi,2);

% Assembly
    asb = Helmholtz(1,1,"gauss",1,domain);
    K = asb.build_stiffness;
    M = asb.build_mass;
    
% Boundary Conditions
    clamp_dofs = domain.extract_boundaries;
    clamp_dofs = clamp_dofs(:,2);
    lm = asb.location_matrix;
    clamp_dofs = lm(:,cell2mat(clamp_dofs(:)));
    clamp_dofs = unique((clamp_dofs(:)));

% Solution
[vecs, omega, solution] = asb.eigensolve(K,M,1,clamp_dofs);
solution = solution{1};
solution.plot_solution(1);
