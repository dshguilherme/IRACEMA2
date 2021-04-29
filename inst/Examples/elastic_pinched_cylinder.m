%% Given Constants
E = 3e6;
vu = 0.3;
P = 1.0;
rho = 300;
L = 600;
t = 3;

%% Geometry

r = rho-t/2;
R = rho+t/2;

point1 = [r 0 0];
point2 = [0 0 r];
center = [0 0 0];
inside_arc = bs_arc(point1, point2, pi/2, center);

point1 = [R 0 0];
point2 = [0 0 R];
outside_arc = bs_arc(point1, point2, pi/2, center);

quarter_annulus = bs_ruled_surface(inside_arc, outside_arc);

domain = bs_extrusion(quarter_annulus, [0 L/2 0]);

%% Assembler
asb = Elastic(E,vu,"gauss",3,domain);
K = asb.build_stiffness;

% Boundary Conditions
boundaries = domain.extract_boundaries;