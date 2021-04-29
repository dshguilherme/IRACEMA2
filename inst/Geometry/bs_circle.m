function geometry_obj = bs_circle(center, radius, plane_vec, plane_point)
assert(length(center) == 3, "ERROR: The center must have 3 coordinates");
assert(length(plane_vec) == 3, "ERROR: Normal vector to the plane must have 3 coordinates");
assert(length(plane_point) == 3, "ERROR: The point of the plane must have 3 coordinates");


in_vec = plane_point - center;
in_vec = in_vec/norm(in_vec);
perp = cross(plane_vec,in_vec);
perp = perp/norm(perp);

B1 = center +radius*in_vec;
B2 = B1 +radius*perp;
B3 = B2 -radius*in_vec;
B4 = B3 -radius*in_vec;
B5 = B4 -radius*perp;
B6 = B5 -radius*perp;
B7 = B6 +radius*in_vec;
B8 = B7 +radius*in_vec;
B9 = B1;

U = [0 0 0 1 1 2 2 3 3 4 4 4]/4;

w = [1 sqrt(2)/2 1 sqrt(2)/2 1 sqrt(2)/2 1 sqrt(2)/2 1];
w = w';
P = [B1;B2;B3;B4;B5;B6;B7;B8;B9];
P = [P,w];
P = num2cell(P,2);

geometry_obj = Geometry(1,{U},P,[2]);

end