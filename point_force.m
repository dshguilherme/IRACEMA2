function f = point_force(x,position,tolerance)
if abs(x-position(1)) < tolerance
    f = [0 -1];
else
    f = [0 0];
end