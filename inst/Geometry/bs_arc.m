function geometry_obj = bs_arc(center, initial_point, theta, normal)

in_vec = initial_point - center;
R = norm(in_vec);

ortho_vec = cross(normal,in_vec);
ortho_vec = ortho_vec/norm(ortho_vec);

P1 = initial_point;
x = R*tan(theta/2);
P2 = initial_point +x*ortho_vec;

mid_vec = P2-center;
ortho_vec = cross(normal,mid_vec);
ortho_vec = ortho_vec/norm(ortho_vec);
P3 = P1 +2*R*sin(theta/2)*ortho_vec;


points = {[P1 1], [P2 cos(theta/2)], [P3 1]};

geometry_obj = Geometry(1,{[0 0 0 1 1 1]},points,[2]);

end
