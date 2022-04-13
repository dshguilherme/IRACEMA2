clearvars
clc
%cantilever beam
E = 1e8; %change later, just testing
vu = 0.33;
L = 10;
alfa = 1;
w = 1;
line = Geometry(1,{[0,0,1,1]},{[0,0,0,1],[L,0,0,1]},1);
line2 = bs_translation(line,[0 w 0]);
domain = bs_ruled_surface(line,line2);
% domain = bs_extrusion(domain,[0 0 w]);

%refinement
Xi = linspace(0,1,36);
Xi = Xi(2:end-1); %no repeating knots
% domain.degree_elevate(1,1);
% domain.degree_elevate(1,2);
%domain.degree_elevate(1,3); %no refine in width
domain.knot_refine(Xi,1);
domain.knot_refine(Xi,2);
%domain.knot_refine(Xi,3);

%asb
asb = Elastic(E,vu,"gauss",2,domain);
K = asb.build_stiffness;
id = asb.id_matrix;
F = zeros(numel(id),1);
P = domain.points;
P = cell2mat(P(:));
f_points = find(P(:,1) == L & P(:,2) == w);
f_dofs = id(f_points,2);
F(f_dofs) = -alfa;



%bcs
x_points = find(P(:,1) == 0);
xyz_dofs = id(x_points,[1 2]);


clamped_dofs = unique(xyz_dofs(:));
g = 0;
[d, F, solution] = asb.dirichlet_linear_solve(K,F,g,clamped_dofs);
solution.plot_solution(2); %y displacement
colormap("jet");