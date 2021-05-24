clearvars
clc
close all

%% Geometry and Material

R = 1;
L = 5;
r = R-0.01/2;
R = R+.01/2;
center = [0 0 0];
line1 = bs_circle(center,r,[0 0 1],[r 0 0]);
line2 = bs_circle(center,R,[0 0 1],[R 0 0]);
cross = bs_ruled_surface(line2,line1);
domain = bs_extrusion(cross,[0 0 L]);

E = 1e3;
vu = 0.33;
pressure = 1.0;

%% Refinement
domain.degree_elevate(1,2);
domain.knot_refine(0.5,2);

domain.degree_elevate(1,3);
domain.degree_elevate([0.1 1/3 2/3 .9],3);

%% Assembly
asb = Elastic(E,vu,"gauss",3,domain);
K = asb.build_stiffness;
F = zeros(length(K),1);

%% Boundary Conditions
id = asb.id_matrix;
bcs = domain.extract_boundaries;
h = @(x) -pressure*[sin(atan(x(2)/x(1))), cos(atan(x(2)/x(1))),0];
F = asb.variable_neumann_bc(F,h,bcs(3,:));
% F = asb.variable_neumann_bc(F,h,bcs(4,:));

clamp_points = [bcs{5,2}; bcs{6,2}];
clamp_dofs = unique(id(:,clamp_points(:)));
g = zeros(numel(clamp_dofs),1);

% Fix points that are coincident
s = size(domain.points,1);
idx1 = 1:s:numel(domain.points);
idx2 = s:s:numel(domain.points);
for i=1:numel(idx1)    
    pdofs = id(:,[idx1(i) idx2(i)]);
    penaltyStiffness = 100000*[1 -1; -1 1];
    px = pdofs(1,:);
    py = pdofs(2,:);
    K(px,px) = K(px,px) +penaltyStiffness;
    K(py,py) = K(py,py) +penaltyStiffness;
end

[d, F, solution] = asb.dirichlet_linear_solve(K,F,g,clamp_dofs);
solution.plot_solution(1)
colormap("jet")
