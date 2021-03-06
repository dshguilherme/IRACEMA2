%% Example 4: Elastic stress of a plate with hole
%% Still under construction!
% Geometric and material
L = 4;
P1 = [-L 0 0 1];
P2 = [-L L 0 1];
P3 = [0 L 0 1];
U = [0 0 0.5 1 1];

line = Geometry(1, {U}, {P3,P2,P1}, [1]);

center = [0 0 0];
initial_point = [0 1 0];
theta = pi/2;
normal = [0 0 1];
arc = bs_arc(center, initial_point, theta, normal);

domain = bs_ruled_surface(line,arc);
domain = bs_extrusion(domain,[0 0 0.01]);

% Refinement
Xi = linspace(0,1,7);
Xi = Xi(2:end-1);
domain.uniform_k_refine(Xi,1);

E = 10e5;
vu = 0.3;

% Assembly

asb = Elastic(E,vu,"gauss",3,domain);
K = asb.build_stiffness;

% Boundary Conditions
boundaries = domain.extract_boundaries;

neumann_boundaries = boundaries(3,:);

theta = @(x,y) atan(-y/x);
r = @(x,y) sqrt(x^2 +y^2);

tau_xx = @(r,theta) 1 - (1/r^2)*(1.5*cos(2*theta) +cos(4*theta)) + ...
                    (1.5/(r^4))*cos(4*theta);
tau_yy = @(r,theta) -(1/r^2)*(0.5*cos(2*theta) -cos(4*theta)) + ...
                    (1.5/(r^4))*cos(4*theta);
tau_xy = @(r,theta) -(1/r^2)*(0.5*sin(2*theta) +sin(4*theta)) + ...
                    (1.5/(r^4))*sin(4*theta);
                
h = @(x) [tau_xx(r(x(1),x(2)),theta(x(1),x(2)));
          tau_yy(r(x(1),x(2)),theta(x(1),x(2)));
          0*tau_xy(r(x(1),x(2)),theta(x(1),x(2)))];

[s1 s1] = size(K);
F = zeros(s1,1);
F = asb.variable_neumann_bc(F,h,neumann_boundaries);
id = asb.id_matrix;


% Dirichlet Boundary Conditions:



dirichlet_boundaries = boundaries([1 2 5 6],:);
x_boundary = dirichlet_boundaries(1,:);
x_points = x_boundary{2};
x_dofs = id(1,x_points);

y_boundary = dirichlet_boundaries(2,:);
y_points = y_boundary{2};
y_dofs = id(2,y_points);

z_boundary = dirichlet_boundaries([3 4],:);
z_points = cell2mat(z_boundary(:,2));
z_dofs = id(3,:);


clamped_dofs = unique([x_dofs(:);y_dofs(:); z_dofs(:)]);
g = zeros(numel(clamped_dofs),1);
[d, F, solution] = asb.dirichlet_linear_solve(K,F,g,clamped_dofs);
[stress, strain] = asb.stress_strain(d);
