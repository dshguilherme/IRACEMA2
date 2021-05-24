function [R, dR, J] = FastShape1D(GeometryObject,IntegrationPoint, ...
    global_basis_index, element_local_mapping, element_ranges, element)

qu = IntegrationPoint(1);
pu = GeometryObject.p(1);
U = GeometryObject.knots{1};

support = global_basis_index(element_local_mapping(:,element),:);

u_range = element_ranges(element,:,1);
u = ((u_range(2)-u_range(1))*qu +(sum(u_range)))/2; % Parent -> Parametric

su = FindSpanLinear(length(U)-pu-2,pu,u,U);

P = GeometryObject.points;
ind = sub2ind(size(P),support(:,1));
ActivePoints = P(ind);
ActivePoints = cell2mat(ActivePoints(:));
Weights = ActivePoints(:,4);
P = ActivePoints(:,1:3);

Basis = DersBasisFun(su,u,pu,1,U);
B = Basis(1,:);
dB = Basis(2,:);
Q = B*Weights;
dQ = dB*Weights;

R = B'.*Weights/Q;
dRdu = Weights.*(Q*dB'-dQ*B')/(Q*Q);
% x = sum((R.*P));
dxdu = sum(P.*dRdu);
dXdU = dxdu';
dUdX = pinv(dXdU);
dudx = dXdU(:,1)';

tmp = element_ranges(element,2,:) - element_ranges(element,1,:);
J = norm(dxdu*tmp(1));
dR = dRdu*dudx;
end