function geometry_obj = bs_circle(center, radius, plane_vec, plane_point)
assert(length(center) == 3, "ERROR: The center must have 3 coordinates");
assert(length(plane_vec) == 3, "ERROR: Normal vector to the plane must have 3 coordinates");
assert(length(plane_point) == 3, "ERROR: The point of the plane must have 3 coordinates");


in_vec = plane_point - center;
in_vec = in_vec/norm(in_vec);
perp = cross(plane_vec,in_vec);
perp = perp/norm(perp);

B1 = center +radius*in_vec;
B3 = center +radius*perp;

P2 = (B3-B1)/2;
v2 = P2-center;
v2 = v2/norm(v2);
B2 = center +radius*v2;

B5 = center -radius*in_vec;

P4 = (B5-B3)/2;
v4 = P4-center;
v4 = v4/norm(v4);
B4 = center +radius*v4;

B7 = center -radius*perp;

P6 = (B7-B5)/2;
v6 = P6-center;
v6 = v6/norm(v6);
B6 = center +radius*v6;

P8 = (B1-B7)/2;
v8 = P8-center;
v8 = v8/norm(v8);
B8 = center +radius*v8;

B9 = B1;

U = [0 0 0 1 1 2 2 3 3 4 4 4]/4;

w = repmat([1 sqrt(2)/2 1],[1,3]);
w = w';
P = [B1;B2;B3;B4;B5;B6;B7;B8;B9];
P = [P,w];
P = num2cell(P,2);

geometry_obj = Geometry(1,{U},P,[2]);

end