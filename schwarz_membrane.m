clearvars
close all
clc

% Geometries 
    L = 1.5;
    h = 1.0;
    o = 0.1;
    b = L*(1-o)/2;
    ell = b+L*o;
    P1 = [0 0 0 1];
    P2 = [ell 0 0 1];
    U = [0 0 1 1];
    
    line1 = Geometry(1, {U}, {P1,P2}, [1]);
    line2 = bs_translation(line1, [0 h 0]);
    omega1 = bs_ruled_surface(line1, line2);
    
    line3 = bs_translation(line1, [b 0 0]);
    line4 = bs_translation(line2, [b 0 0]);
    omega2 = bs_ruled_surface(line3, line4);
    
% Refinement
    Xi = linspace(0,1,15);
    Xi = Xi(2:end-1);
    omega1.degree_elevate(1,1);
    omega2.degree_elevate(1,1);
    omega1.degree_elevate(1,2);
    omega2.degree_elevate(1,2);
    omega1.knot_refine(Xi,1);
    omega2.knot_refine(Xi,1);
    omega1.knot_refine(Xi,2);
    omega2.knot_refine(Xi,2);

% Assemblies
    asb1 = Helmholtz(1,1,"gauss",1,omega1);
    [K1, M1] = asb1.build_matrices;
    
    asb2 = Helmholtz(1,1,"gauss",1, omega2);
    [K2, M2] = asb2.build_matrices;

% Dirichlet Boundary Conditions
    boundaries = omega1.extract_boundaries;
    boundary_1 = boundaries([1 3 4],:);
    clamp_1_elements = boundary_1(:,2);
    lm_1 = asb1.location_matrix;
    clamp_1_dofs = lm_1(:,cell2mat(clamp_1_elements(:)));
    clamp_1_dofs = unique((clamp_1_dofs(:)));
    
    boundaries = omega2.extract_boundaries;
    boundary_2 = boundaries([2 3 4],:);
    clamp_2_elements = boundary_2(:,2);
    lm_2 = asb2.location_matrix;
    clamp_2_dofs = lm_2(:,cell2mat(clamp_2_elements(:)));
    clamp_2_dofs = unique((clamp_2_dofs(:)));
    
% Initial Solutions
    [u_1, w_1, sol_1] = asb1.eigensolve(K1,M1,1,clamp_1_dofs);
    sol_1 = sol_1{1};
    figure(1)
    sol_1.plot_solution(1);
    title('Left Domain Initial Solution');    
    [u_2, w_2, sol_2] = asb2.eigensolve(K2,M2,1,clamp_2_dofs);
    sol_2 = sol_2{1};
    figure(2)
    sol_2.plot_solution(1);
    title('Right Domain Initial Solution');

% Building inverse mappings
    omega1_inverse_mapping = @(x) [x(1)/ell, x(2)/h];
    omega2_inverse_mapping = @(x) [(x(1)-b)/ell, x(2)/h];
    solution = sol_2;

    C = 5e3;
    gamma = 1;
    KK1 = K1;
    for i=1:1
% First iterative step: weak enforcement of dirichlet BCs on Omega 1
    g1 = @(x) solution.eval_solution_value(omega2_inverse_mapping(x));

    boundaries1 = omega1.extract_boundaries;
    bdries = boundaries1(2,:);
    [KK, FF] = asb1.weak_dirichlet(g1,bdries,C,gamma);
    K1 = K1+KK;
    F1 = FF;
    [d, F1, solution] = asb1.dirichlet_linear_solve(K1,F1,0,clamp_1_dofs);
    sol1 = solution;
% Second iterative step: weak enforcement of dirichlet BCs on Omega 1
    g2 = @(x) solution.eval_solution_value(omega1_inverse_mapping(x));
    C = 1e6;
    gamma = 1;
    boundaries2 = omega2.extract_boundaries;
    bdries = boundaries(1,:);
    [KK, FF] = asb2.weak_dirichlet(g2,bdries,C,gamma);
    K2 = K2+KK;
    F2 = FF;
    [d, F2, solution] = asb2.dirichlet_linear_solve(K2,F2,0,clamp_2_dofs);
    sol2 = solution;
    end
    figure(3)
        sol1.plot_solution(1)
        hold on
        sol2.plot_solution(1)

