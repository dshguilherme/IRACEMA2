function [nn, nU, nP] = SurfaceKnotInsert(n,p,U,P,xi,dir)
assert(length(xi) == 1, "ERROR: Must insert only one knot");
nu = n-1;
pu = p;
k = FindSpanLinear(nu,pu,xi,U);
k = k+1;
s = size(P);
s(dir) = s(dir)+1;

notdir = setdiff([1 2],dir);
idx1 = cell(2,1);
idx1{dir} = 1:k-pu;
idx1{notdir} = 1:s(notdir);

idx2 = cell(2,1);
idx2{dir} = k+1:s(dir);
idx2{notdir} = 1:s(notdir);

idx3 = cell(2,1);
idx3{notdir} = 1:s(notdir);

Qw = cell(size(P));
s2 = size(P);
Pw = P(:);
Pw = cell2mat(Pw);
Pw(:,1:3) = Pw(:,1:3).*Pw(:,4);
Pw = num2cell(Pw,2);
Pw = reshape(Pw,s2);

Qw(idx1{1},idx1{2}) = Pw(idx1{1},idx1{2});
idxx = idx2;
idxx{dir} = idxx{dir}-1;
Qw(idx2{1},idx2{2}) = Pw(idxx{1},idxx{2});

for i = k-pu+1:k
    alpha = (xi-U(i))/(U(i+pu)-U(i));
    idx3{dir} = i;
    Qw1 = Pw(idx3{1},idx3{2});
    idx3{dir} = i-1;
    Qw2 = Pw(idx3{1},idx3{2});
    Qw3 = Qw1;
    for j=1:length(Qw1)
        Qw3{j} = alpha*Qw1{j} +(1-alpha)*Qw2{j};
    end
    idx3{dir} = i;
    Qw(idx3{1},idx3{2}) = Qw3;
end
nU = [U(1:k) xi U(k+1:end)];
nP = cell2mat(Qw(:));
nP(:,1:3) = nP(:,1:3)./nP(:,4);
nP = num2cell(nP,2);
nP = reshape(nP,size(Qw));
nn = length(nU)-pu-1;


end