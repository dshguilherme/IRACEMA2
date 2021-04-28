function geometry_obj = bs_line(point1, point2)
assert(length(point1) == 3), "ERROR: point1 input size must be 3");
assert(length(point2) == 3), "ERROR: point2 input size must be 3");
U = [0 0 1 1];
P1 = [point1 1];
P2 = [point2 1];

geometry_obj = Geometry(1, {U}, {P1, P2}, [1]);
end