function geometry_obj = bs_arc(point1, point2, theta, plane_point)
assert(length(point1) == 3, "ERROR: point1 input size must be 3");
assert(length(point2) == 3, "ERROR: point2 input size must be 3");
assert(theta <= pi, "ERROR: theta must be between 0 and pi");
assert(norm(plane_point -point1) > 0, "ERROR: point in plane must not be coincident with p1 or p2");
assert(norm(plane_point -point2) > 0, "ERROR: point in plane must not be coincident with p1 or p2");

B1 = point1;
B2 = point2;
B3 = plane_point;

P1 = (B2-B1)/2;

in_vec1 = B2-P1;
in_vec2 = B3-P1;

N = cross(in_vec2,in_vec1);
n = N/norm(N);

perpendicular = cross(in_vec1,n);
perpendicular = perpendicular/norm(perpendicular);

ell = norm(B2-B1);

B4 = P1 + (ell/2)*tan(theta/2)*perpendicular;

points = {[B1 1], [B4 cos(theta/2)], [B2 1]};

geometry_obj = Geometry(1,{[0 0 0 1 1 1]},points,[2]);

end