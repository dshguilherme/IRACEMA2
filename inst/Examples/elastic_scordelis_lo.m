%% Given Constants
E = 4.32e8;
vu = 0.0;
g = [0 0 -90];
R = 25;
L = 50;
phi = 40*pi/180;
t = 0.25;

%% Geometry
theta = pi/2 - phi;
point1 = [0 0 R+t/2];
point2 = [(R+t/2)*cos(theta) 0 (R+t/2)*sin(theta)];

upper = bs_arc(point1, point2, phi, [0 0 0]);

point1 = [0 0 R-t/2];
point2 = [(R-t/2)*cos(theta) 0 (R-t/2)*sin(theta)];

lower = bs_arc(point1,point2, phi, [0 0 0]);

cross_section = bs_ruled_surface(upper,lower);

domain = bs_extrusion(cross_section,[0 25 0] );

%% Assembler
asb = Elastic(E,vu,"gauss",3,domain);
K = asb.build_stiffness;
F = asb.constant_force(g);

% Boundary Conditions
boundaries = domain.extract_boundaries;

