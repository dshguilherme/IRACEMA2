clearvars
clc
close all
%% Geometry and Material

P =[
    2.0000         0         0    1.0000;
    4.0000         0         0    0.7071;
    4.0000    2.0000         0    1.0000;
         0         0         0    0.7071;
    2.0000    2.0000         0    1.0000;
    4.0000    4.0000         0    0.7071;
         0    2.0000         0    1.0000;
         0    4.0000         0    0.7071;
    2.0000    4.0000         0    1.0000;
    2.0000         0    0.0200    1.0000;
    4.0000         0    0.0200    0.7071;
    4.0000    2.0000    0.0200    1.0000;
         0         0    0.0200    0.7071;
    2.0000    2.0000    0.0200    1.0000;
    4.0000    4.0000    0.0200    0.7071;
         0    2.0000    0.0200    1.0000;
         0    4.0000    0.0200    0.7071;
    2.0000    4.0000    0.0200    1.0000];
P = num2cell(P,2);
P = reshape(P,[3 3 2]);
p = [2 2 1];
U = [0 0 0 1 1 1];
V = [0 0 0 1 1 1];
W = [0 0 1 1];

domain = Geometry(3,{U,V,W},P,p);
E = 30e6;
vu = 0.2;
rho = 2320;

% Refinement
domain.degree_elevate(2,1);
domain.degree_elevate(3,2);
domain.degree_elevate(1,3);

% Assembly
asb = Elastodynamic(rho,E,vu,"gauss",3,domain);
[K,M] = asb.build_matrices;

% Boundary Conditions
id = asb.id_matrix;
bcs = domain.extract_boundaries;
clamp_points = cell2mat(bcs(1:4,2));
clamp_points = unique(clamp_points);
clamp_dofs = id(:,clamp_points);
clamp_dofs = clamp_dofs(:);

% Solve
[vecs, omega, solution_cell] = asb.eigensolve(K,M,10,clamp_dofs);



