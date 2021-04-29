function geometry_obj = bs_rotation(geometry, axis, theta)
RM = [cos(theta) + axis(1)*axis(1)*(1-cos(theta)), axis(1)*axis(2)*(1-cos(theta)) - axis(3)*sin(theta), axis(1)*axis(3)*(1-cos(theta)) + axis(2)*sin(theta);
                   axis(2)*axis(1)*(1-cos(theta)) + axis(3)*sin(theta), cos(theta) + axis(2)*axis(2)*(1-cos(theta)), axis(2)*axis(3)*(1-cos(theta)) - axis(1)*sin(theta);
                   axis(3)*axis(1)*(1-cos(theta)) - axis(2)*sin(theta), axis(3)*axis(2)*(1-cos(theta)) + axis(1)*sin(theta), cos(theta) + axis(3)*axis(3)*(1-cos(theta))];

P = geometry.points;
s = size(P);
B = cell(s);
P = cell2mat(P(:));
for i=1:prod(s)
    B{i} = [rot*P(i,1:3) P(i,4)];
end

geometry_obj = Geometry(geometry.rank, geometry.knots, B, geometry.p);

end