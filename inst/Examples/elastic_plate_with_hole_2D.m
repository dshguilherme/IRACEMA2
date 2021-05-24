clearvars
clc
%% Example 4: 2D Plate with circular hole

% Geometric and Material
L = 4;
P1 = [-L 0 0 1];
P2 = [-L L 0 1];
P3 = [0 L 0 1];
U = [0 0 0.5 1 1];

line = Geometry(1, {U}, {P3,P2,P1}, [1]);
R = 1;
center = [0 0 0];
initial_point = [0 R 0];
theta = pi/2;
normal = [0 0 1];
arc = bs_arc(center, initial_point, theta, normal);

domain = bs_ruled_surface(line,arc);

% Refinement
Xi = linspace(0,1,7);
Xi = Xi(2:end-1);
domain.uniform_k_refine(Xi,3);

E = 1000;
vu = 0.3;

% Assembly

asb = Elastic(E,vu,"gauss",2,domain);
K = asb.build_stiffness;

% Boundary Conditions
boundaries = domain.extract_boundaries;

x_sym_bc = boundaries(1,:);
y_sym_bc = boundaries(2,:);
stress_bc = boundaries(3,:);
zero_stress_bc = boundaries(4,:);

% Neumann BCs
id = asb.id_matrix;
[s1 s1] = size(K);
F = zeros(s1,1);
% Top domain
fasb = Assembler("gauss",2,stress_bc{1});
[qp, qw] = fasb.quad_rule;
[global_basis_index, element_local_mapping, element_ranges] = ...
    GetConnectivityArrays(stress_bc{1});
[~, lm] = BuildGlobalLocalMatrices(element_local_mapping,2);

n_quad = length(qw);
ndof = max(max(element_local_mapping))*2;
[nel_dof, nel] = size(element_local_mapping);
Fi = zeros(ndof,1);
for e=1:nel
    F_e = zeros(2,nel_dof);
    for n=1:n_quad
        q = qp(n,:);
        [R, ~, J] = FastShape(stress_bc{1}, q, global_basis_index, ...
            element_local_mapping, element_ranges, e);
        Jmod = abs(J*qw(n));
        u = 0.5*q +0.5;
        x = stress_bc{1}.eval_point(u);
        st = exact_plate_hole(x,1);
        
        if x(1) > -4
            tx = st(3);
            ty = st(2);
        else
            tx = -st(1);
            ty = -st(3);
        end
        F_e(1,:) = F_e(1,:) +R'*tx*Jmod;
        F_e(2,:) = F_e(2,:) +R'*ty*Jmod;
    end
    idx = lm(:,e)';
    Fi(idx) = Fi(idx) +F_e(:);
end
fpoints = stress_bc{2};
fdofs = id(:,fpoints);
F(fdofs(:)) = F(fdofs(:)) + Fi;
    

% Dirichlet BCs

x_points = x_sym_bc{2};
x_dofs = id(1,x_points);

y_points = y_sym_bc{2};
y_dofs = id(2,y_points);

clamped_dofs = unique([x_dofs(:); y_dofs(:)]);

g = zeros(numel(clamped_dofs,1));

% Penalty to enforce the overlapping control points
% to have the same displacement
P = domain.points;
P = cell2mat(P(:));
[s1, ~] = size(domain.points);

ppoint = find(P(:,1) == -L & P(:,2) == L);
for i=ppoint(1):s1:length(P)
    pdofs = id(:,[i i+1]);
    penaltyStiffness = 1000*[1 -1; -1 1];
    px = pdofs(1,:);
    py = pdofs(2,:);
    K(px,px) = K(px,px) +penaltyStiffness;
    K(py,py) = K(py,py) +penaltyStiffness;
end

[d, F, solution] = asb.dirichlet_linear_solve(K,F,g,clamped_dofs);
solution.plot_solution(1)
[stress sd] = asb.calculate_stresses(d);
% figure
% clf
% domain.plot_geo
asb2 = Assembler("gauss",length(asb.D),domain);
solution = Solution(asb2,sd);
solution.plot_solution(1)
% colormap("jet")
