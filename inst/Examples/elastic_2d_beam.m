clearvars
clc
%% Example: 2D Beam
% Material Constants
rho = 1; % 
E = 1e3;
vu = 0.33;
omega = 0;
alpha = 1;
beta = 1e-8;
% Geometry
L = 20;
b = 1;

P1 = [0 0 0];
P2 = [L 0 0];

line = bs_line(P1,P2);
line2 = bs_translation(line,[0 b 0]);

domain = bs_ruled_surface(line,line2);

% Refinement
domain.degree_elevate(1,1);
domain.degree_elevate(1,2);

domain.knot_refine([1/4 0.5 3/4],2);

Xi = linspace(0,1,101);
Xi = Xi(2:end-1);
domain.knot_refine(Xi,1);

% Assembly
asb = Elastodynamic(rho,E,vu,"gauss",2,domain);
K = asb.build_stiffness;
M = asb.build_mass;

% % BCs
id = asb.id_matrix;
boundaries = domain.extract_boundaries;
clamp_dofs = boundaries(1,2);
lm = asb.location_matrix;
clamp_dofs = lm(:,cell2mat(clamp_dofs(:)));
clamp_dofs = unique((clamp_dofs(:)));

% Damping
C = alpha*M +beta*K;
Kd = @(w) K +1i*C -(w*w)*M;

% Solution
F = zeros(length(K),1);

f = @(x) point_force(L,x,0.05);
F = asb.force_vector(f);
g = 0;

valores = zeros(101,1);
for i=1:101
[d, F, solution] = asb.dirichlet_linear_solve(Kd(i/5-1/5),F,g,clamp_dofs);
uu = solution.eval_solution_value([1 0.5]);
valores(i) = abs(uu(2));
end
% ylim([-10 10]);
semilogx(0:1/5:20,valores)
