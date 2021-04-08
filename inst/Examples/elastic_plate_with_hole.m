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
arc = geo_arc(center,1,'xy',2);

domain = geo_ruled(line,arc);
domain = geo_extrusion(domain,[0 0 0.01]);

E = 10e5;
vu = 0.3;

% Assembly

asb = Elastic(E,vu,"gauss",2,domain);
K = asb.build_stiffness;

% Boundary Conditions
boundaries = domain.extract_boundaries;

theta = @(x,y) atan(-y/x);
r = @(x,y) sqrt(x^2 +y^2);

tau_xx = @(r,theta) 1 - (1/r^2)*(1.5*cos(2*theta) +cos(4*theta)) + ...
                    (1.5/(r^4))*cos(4*theta);
tau_yy = @(r,theta) -(1/r^2)*(0.5*cos(2*theta) -cos(4*theta)) + ...
                    (1.5/(r^4))*cos(4*theta);
tau_xy = @(r,theta) -(1/r^2)*(0.5*sin(2*theta) +sin(4*theta) + ...
                    (1.5/(r^4))*sin(4*theta);


% Dirichlet Boundary Conditions: