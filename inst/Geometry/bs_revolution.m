function geometry_obj = bs_revolution(geometry, axis, theta)
if theta <= pi/2
        narcs = 1;
        W = [0 0 0 1 1 1];
    elseif theta <= pi
        narcs = 2;
        W = [0 0 0 0.5 0.5 1 1 1];
    elseif theta <= 3*pi/2
        narcs = 3;
        W = [0 0 0 1 1 2 2 3 3 3]/3;
    elseif theta <= 2*pi
        narcs = 4;
        W = [0 0 0 1 1 2 2 3 3 4 4 4]/4;
    else
        error("ERROR: theta must be between 0 and 2*pi");
end

P = geometry.points;
if geometry.rank == 1
    s = [length(P), 2*narcs+1];
else
    s = [size(P), 2*narcs+1];
end
P = cell2mat(P(:));
theta_vec = linspace(0,theta,2*narcs+1);
dtheta = theta/narcs;

B = cell(s);
B = B(:);
for i=1:2*narcs+1
    rotation = bs_rotation(geometry, axis, theta_vec(i));
    P = rotation.points;
    P = cell2mat(P(:));
    P(:,4) = P(:,4)*weights(i);
    P = num2cell(P,2);
    B(1+(i-1)*numel(P):numel(P)+(i-1)*numel(P)) = P;
end
B = reshape(B,s);

geometry_obj = Geometry(geometry.rank+1,[geometry.knots {W}],B,[geometry.p 2]);    
end