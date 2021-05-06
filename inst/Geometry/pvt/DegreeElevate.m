function [nn, newU, Q] = DegreeElevate(n,p,U,P)
adds = unique(U);
newU = sort([U adds]);

nn = length(newU)-(p+1)-1;

Pw = cell2mat(P(:));
Pw(:,1:3) = Pw(:,1:3).*Pw(:,4);
Px = Pw(:,1);
Py = Pw(:,2);
Pz = Pw(:,3);
Pww = Pw(:,4);

samples = linspace(0,1,nn);
NP = zeros(nn,n);
NP(1) = 1;
NP(end) = 1;
MP = zeros(nn,nn);
MP(1) = 1;
MP(end) = 1;
for i=2:nn-1
    u = samples(i);
    s = FindSpanLinear(n-1,p,u,U);
    N = DersBasisFun(s,u,p,0,U);
    NP(i,s-p+1:s+1) = N;
    su = FindSpanLinear(nn-1,p+1,u,newU);
    M = DersBasisFun(su,u,p+1,0,newU);
    MP(i,su-(p+1)+1:su+1) = M;
end
T = MP\NP;
Qw = [T*Px, T*Py, T*Pz, T*Pww];
Qw(:,1:3) = Qw(:,1:3)./Qw(:,4);
Q = num2cell(Qw,2);
Q = Q';

end