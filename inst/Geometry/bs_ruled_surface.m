function geometry_obj = bs_ruled_surface(curve1, curve2)
p1 = curve1.p;
p2 = curve2.p;

p = max(p1, p2);

if p1 > p2
    curve2.degree_elevate(p1-p2,1);
elseif p2 > p1
    curve1.degree_elevate(p2-p1,1);
end

nu = max(curve1.n,curve2.n);

U1 = curve1.knots{1};
U2 = curve2.knots{1};

tmp2 = ~ismember(U1,U2);
knots_to_insert2 = U1(tmp2);

tmp1 = ~ismember(U2,U1);
knots_to_insert1 = U2(tmp1);

if knots_to_insert2(1) ~= 0
    curve2.knot_refine(knots_to_insert2,1);
end
if knots_to_insert1(1) ~= 0
    curve1.knot_refine(knots_to_insert1,1);
end

B1 = curve1.points;
B2 = curve2.points;

U = U1;
V = [0 0 1 1];

S = cell(curve1.n(1) +1,2);
for i = 1:curve1.n(1) +1
    S{i,1} = B1{i};
    S{i,2} = B2{i};
end

geometry_obj = Geometry(2, {U,V}, S, [p, 1]);

end