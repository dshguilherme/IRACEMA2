function geometry_obj = bs_scaling(geometry, vector)
assert(length(vector) == 3, "ERROR: scaling vector must be tri-dimensional");

S = diag(vector);
P = geometry.points;
s = size(P);
B = cell(s);
P = cell2mat(P(:));
for i=1:prod(s)
    B{i} = [S*P(i,1:3) P(i,4)];
end

geometry_obj = Geometry(geometry.rank, geometry.knots, B, geometry.p);

end