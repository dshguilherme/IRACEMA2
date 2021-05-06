function [nn, newU, Q] = SurfaceDegreeElevate(n,p,U,P,dir)
adds = unique(U);
newU = sort([U adds]);

nn = length(newU)-(p+1)-1;

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
[s1, s2] = size(P);
switch dir
    case 1
        for i=1:s2
            AP = P(:,i);
            APw = cell2mat(AP(:));
            APw(:,1:3) = APw(:,1:3).*APw(:,4);
            Px = APw(:,1);
            Py = APw(:,2);
            Pz = APw(:,3);
            Pww = APw(:,4);
            Qw = [T*Px, T*Py, T*Pz, T*Pww];
            Qw(:,1:3) = Qw(:,1:3)./Qw(:,4);
            Qw = num2cell(Qw,2);
            Q(:,i) = Qw;
        end
    case 2
        for i=1:s1
           AP = P(i,:);
            APw = cell2mat(AP(:));
            APw(:,1:3) = APw(:,1:3).*APw(:,4);
            Px = APw(:,1);
            Py = APw(:,2);
            Pz = APw(:,3);
            Pww = APw(:,4);
            Qw = [T*Px, T*Py, T*Pz, T*Pww];
            Qw(:,1:3) = Qw(:,1:3)./Qw(:,4);
            Qw = num2cell(Qw,2);
            Q(i,:) = Qw;
        end
end

end