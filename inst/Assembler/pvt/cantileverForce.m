function f = cantileverForce(x)
fx = 0;
if x(1) > 0.9 && (x(2) < 0.75 || x(2) > 0.5) 
    fy = -1;
else
    fy = 0;
end
f = [fx fy 0];
end