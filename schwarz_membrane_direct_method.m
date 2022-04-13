clearvars
close all
clc

%% Geometries 
    L = 1.5;
    h = 1.0;
    o = 0.1;
    b = L*(1-o)/2;
    ell = b+L*o;
    P1 = [0 0 0 1];
    P2 = [ell 0 0 1];
    P3 = [L 0 0 1];
    U = [0 0 1 1];
    
    line1 = Geometry(1, {U}, {P1,P2}, [1]);
    line2 = bs_translation(line1, [0 h 0]);
    omega1 = bs_ruled_surface(line1, line2);
    
    line3 = bs_translation(line1, [b 0 0]);
    line4 = bs_translation(line2, [b 0 0]);
    omega2 = bs_ruled_surface(line3, line4);
    
    line5 = Geometry(1, {U}, {P1,P3}, [1]);
    line6 = bs_translation(line5, [0 h 0]);
    omega = bs_ruled_surface(line5, line6);
    
%% Refinement
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
    
% Coarser Omega
    Xi = linspace(0,1,10);
    Xi = Xi(2:end-1);
    omega.degree_elevate(1,1);
    omega.degree_elevate(1,2);
    omega.knot_refine(Xi,1);
    omega.knot_refine(Xi,2);
    
%% Assemblies
    asb1 = Helmholtz(1,1,"gauss",1,omega1);
    [K1, M1] = asb1.build_matrices;
    
    asb2 = Helmholtz(1,1,"gauss",1, omega2);
    [K2, M2] = asb2.build_matrices;

    
    asb = Helmholtz(1,1,"gauss",1,omega);
    [K, M] = asb.build_matrices;
    
%% Dirichlet Boundary Conditions
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
    
    boundaries = omega.extract_boundaries;
    clamp_elements = boundaries(:,2);
    lm = asb.location_matrix;
    clamp_dofs = lm(:,cell2mat(clamp_elements(:)));
    clamp_dofs = unique((clamp_dofs(:)));
    
%% Initial Solutions
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
    
    [u, w, sol] = asb.eigensolve(K,M,1,clamp_dofs);
    sol = sol{1};
    figure(3)
    sol.plot_solution(1)
    title('Whole Domain Initial Solution (Coarser)');

%% Mappings
    omega1_inverse_mapping = @(x) [x(1)/ell, x(2)/h];
    omega2_inverse_mapping = @(x) [(x(1)-b)/ell, x(2)/h];
    omega_inverse_mapping = @(x) [x(1)/L, x(2)/h];
   
%% Schwarz Algorithm
    C = 1e4;
    gamma = 1;
    boundaries1 = omega1.extract_boundaries;
    boundaries2 = omega2.extract_boundaries;
    solution2 = sol_2;
    % Step 1: Project the solution of Omega 2 at Omega 1
        g1 = @(x) solution2.eval_solution_value(omega2_inverse_mapping(x));
    % Step 2: Project the derivative of the Omega 2 solution at Omega 1
    % Step 3: Build and solve the new problem
        [KK1, FF1] = asb1.weak_dirichlet(g1, boundaries1(2,:), C,gamma);
        KK1 = K1+KK1;
        [d, FF1, solution1] = asb1.helmholtz_linear_solve(KK1,M1,FF1,w, ...
                                                          0, clamp_1_dofs);
       figure(4)
       solution1.plot_solution(1)
    % Step 4: Project the solution of Omega 1 i+1 at Omega 2 i
        g2 = @(x) solution1.eval_solution_value(omega1_inverse_mapping(x));
    % Step 5: Project the derivative of the Omega 1 i+1 solution at Omega 2
        [KK2, FF2] = asb2.weak_dirichlet(g2, boundaries2(1,:), C, gamma);
        KK2 = K2+KK2;
        [d, FF2, solution2] = asb2.helmholtz_linear_solve(KK2,M2,FF2,w, ...
                                                           0,clamp_2_dofs);
        figure(5)
        solution2.plot_solution(1)
    % Step 6: Check for convergence criteria