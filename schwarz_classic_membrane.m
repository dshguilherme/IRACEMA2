clearvars
close all
clc

% Geometries
center = [0,0,0];
radius = 1;
plane_point = [0, 1, 0];
plane_vec = [0, 0, 1];

circle = bs_circle(center, radius, plane_vec, plane_point);

P1 = [0 radius 0];
P2 = [radius 0 0];
line1 = bs_line(P1, P2);

P3 = [2*radius radius 0];
P4 = [2*radius 0 0];
line2 = bs_line(P3, P4);

rectangle = bs_ruled_surface(line1, line2);
sing_line = bs_line([0, 0, 0], [0, 0, 0]);
circle_2 = bs_ruled_surface(circle, sing_line);

P = cell(1,13);
P{1} = [0 1 0 1];
P{2} = [2 1 0 1];
P{3} = [2 1 0 1];
P{4} = [2 0 0 1];
P{5} = [2 0 0 1];
P{6} = [1 0 0 1];
P{7} = [1 0 0 1];
P{8} = [1 -1 0 1/sqrt(2)];
P{9} = [0 -1 0 1];
P{10} = [-1 -1 0 1/sqrt(2)];
P{11} = [-1 0 0 1];
P{12} = [-1 1 0 1/sqrt(2)];
P{13} = P{1};
U = [0 0 0 1 1 2 2 3 3 4 4 5 5 6 6 6];
geo = Geometry(1,{U},P,2);

domain = bs_ruled_surface(geo,line1);

omega = domain;
omega1 = circle_2;
omega2 = rectangle;

%% Refinement
    Xi = linspace(0,1,35);
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
    boundary_1 = boundaries(3,:);
    elm = omega1.element_local_mapping;
    elements = 1:size(elm,2);
    clamp_points = boundary_1(:,4);
    clamp_points = cell2mat(clamp_points(:));
    clamp_1_elements = elements(any(ismember(elm,clamp_points,"legacy")))';
    lm_1 = asb1.location_matrix;
    clamp_1_dofs = lm_1(:,(clamp_1_elements(:)));
    clamp_1_dofs = unique((clamp_1_dofs(:)));
    
    boundaries = omega2.extract_boundaries;
    boundary_2 = boundaries([1 2 4],:);
    clamp_2_elements = boundary_2(:,2);
    lm_2 = asb2.location_matrix;
    clamp_2_dofs = lm_2(:,cell2mat(clamp_2_elements(:)));
    clamp_2_dofs = unique((clamp_2_dofs(:)));
    
    boundaries = omega.extract_boundaries;
    clamp_elements = boundaries(3,2);
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


    