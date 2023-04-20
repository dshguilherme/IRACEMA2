function domain = bs_rectangle(L,h)
P1 = [0 0 0 1];
P2 = [L 0 0 1];
U = [0 0 1 1];
line1 = Geometry(1, {U}, {P1, P2}, [1]);
line2 = bs_translation(line1,[0 h 0]);
domain = bs_ruled_surface(line1, line2);
end