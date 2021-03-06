function [nn nU, nP] = KnotInsert(nu,pu,U,P,xi)
assert(length(xi) == 1,"ERROR: Must insert only one knot")
nu = nu-1;
k = FindSpanLinear(nu,pu,xi,U);
k = k+1;
Pw = P(:);
Pw = cell2mat(Pw);
Pw(:,1:3) = Pw(:,1:3).*Pw(:,4);
[s1, s2] = size(Pw);
Qw = zeros(s1+1,s2);
Qw(1:k-pu,:) = Pw(1:k-pu,:);
Qw(k+1:end,:) = Pw(k:end,:);
for i = k-pu+1:k
    alpha = (xi-U(i))/(U(i+pu)-U(i));
    Qw(i,:) = alpha*Pw(i,:) +(1-alpha)*Pw(i-1,:);
end
nU = [U(1:k) xi U(k+1:end)];
Qw(:,1:3) = Qw(:,1:3)./Qw(:,4);
nP = num2cell(Qw,2);
nP = reshape(nP,[1,length(P)+1]);
nn = length(nU)-pu-1;
end