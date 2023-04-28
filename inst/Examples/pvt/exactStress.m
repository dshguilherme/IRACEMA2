function sigma = exactStress(Tx,R,x)

[theta, r] = cart2pol(x(1),x(2));

T = 0.5*Tx;
rat = R/r;

s_rr = T*(1-rat^2) +T*(1 -4*rat^2 +3*rat^4)*cos(2*theta);
s_tt = T*(1+rat^2) -T*(1 +3*rat^4)*cos(2*theta);
s_rt = -T*(1 +2*rat^2 -3*rat^4)*sin(2*theta);

s_xx = s_rr*(cos(theta)^2) +s_tt*(sin(theta)^2) -s_rt*sin(2*theta);
s_yy = s_rr*(sin(theta)^2) +s_tt*(cos(theta)^2) +s_rt*sin(2*theta);
s_xy = sin(theta)*cos(theta)*(s_rr -s_tt) +s_rt*cos(2*theta);

sigma = [s_xx; s_yy; s_xy];

end