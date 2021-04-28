function geometry_obj = bs_extrusion(plane, vector)
assert(plane.rank == 2, "ERROR: You must select a plane object to extrude")
S = plane.points;
U = plane.knots{1};
V = plane.knots{2};
pu = plane.p(1);
pv = plane.p(2);
nu = plane.n(1);
nv = plane.n(2);

vec = [vector 0];
B = S;
for i=1:nu+1
    for j=1:nv+1
        B{i,j,2} = B{i,j,1} +vec;
    end
end
geometry_obj = Geometry(3,{U,V,[0 0 1 1]},B,[pu pv 1]);
end