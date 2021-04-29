function geometry_obj = bs_extrusion(surface, vector)
assert(surface.rank == 2, "ERROR: You must select a surface object to extrude")
S = surface.points;
U = surface.knots{1};
V = surface.knots{2};
pu = surface.p(1);
pv = surface.p(2);
nu = surface.n(1);
nv = surface.n(2);

vec = [vector 0];
B = S;
for i=1:nu+1
    for j=1:nv+1
        B{i,j,2} = B{i,j,1} +vec;
    end
end
geometry_obj = Geometry(3,{U,V,[0 0 1 1]},B,[pu pv 1]);
end