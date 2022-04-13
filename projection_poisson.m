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
    Xi = linspace(0,1,20);
    Xi = Xi(2:end-1);
    domain.degree_elevate(1,1);
    domain.degree_elevate(1,2);
    domain.knot_refine(Xi,1);
    domain.knot_refine(Xi,2);

% Assembly
    asb = Poisson(1,"gauss",1,domain);

% 1D L2_projector
    Pi = asb.L2_projector(1);
   
% Function to be projected
    f = @(x) sin(pi*x(1)/1.5)*sin(pi*x(2));
    F = asb.project_function(f,1);

% Dirichlet Boundary Conditions
    clamp_dofs = domain.extract_boundaries;
    clamp_dofs = clamp_dofs(:,2);
    lm = asb.location_matrix;
    clamp_dofs = lm(:,cell2mat(clamp_dofs(:)));
    clamp_dofs = unique((clamp_dofs(:)));

% Solve
    d = zeros(size(F));
%     d = Pi\F;
    free_dofs = setdiff(1:length(d),clamp_dofs);
    d(clamp_dofs) = 0;
    F(free_dofs) = F(free_dofs) - Pi(free_dofs,clamp_dofs)*d(clamp_dofs);
    d(free_dofs) = Pi(free_dofs,free_dofs)\F(free_dofs);
    solution = Solution(asb, d);
    solution.asb = asb;

% Post-process
    solution.plot_solution(1);
 
    