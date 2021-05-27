function s = plate_stress(x,Tx)
r = sqrt(x(1)^2 +x(2)^2);
t = atan(x(2)/x(1));

s_rr = (0.5*Tx)*(1-1/r^2) + (0.5*Tx)*(1 -4/r^2 +3/r^4)*cos(2*t);
s_tt = (0.5*Tx)*(1+1/r^2) - (0.5*Tx)*(1+3/r^4)*cos(2*t);
s_rt = -(0.5*Tx)*(1 +2/r^2 -3/r^4)*sin(2*t);

c = cos(t);
s = sin(t);

s_xy = ((c^2 -s^2)*s_rt +c*s*(s_rr -s_tt))/(c^2 -2*c*c*s*s +s^4 -4*c*s);
s_yy = 0.5*(s_rr +s_tt -((c*c -s*s)*(s_xy - s_rt)/(c*s)));
s_xx = s_rr +s_tt -s_yy;


s = [s_xx s_yy s_xy];


end