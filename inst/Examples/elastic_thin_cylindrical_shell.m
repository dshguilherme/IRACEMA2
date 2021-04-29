%% Given Constants
rho = 1;
t = 0.01;
L = 5;
E = 1e6;
vu = 0.3;

%% Geometry

r = rho-t/2;
R = rho+t/2;
center = [0 0 0];
normal = [0 0 1];
point_in_plane = [1 0 0];

inner_circle = bs_circle(center, r, normal, point_in_plane);
outer_circle = bs_circle(center, R, normal, point_in_plane);

annulus = bs_ruled_surface(inner_circle, outer_circle);

domain = bs_extrusion(annulus, [0 0 L]);

%% Assembler
asb = Elastic(E,vu,"gauss",3,domain);
K = asb.build_stiffness;


