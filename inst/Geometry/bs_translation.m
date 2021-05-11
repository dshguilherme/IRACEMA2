function geometry_obj = bs_translation(geometry, vector)
assert(length(vector) == 3, "ERROR: input vector must be tridimensional");

vector = [vector 0];
P = geometry.points;
s = size(P);
B = cell(s);
P = cell2mat(P(:));
for i=1:prod(s)
    B{i} = P(i,:) +vector;
end
geometry_obj = Geometry(geometry.rank, geometry.knots, B, geometry.p);

end