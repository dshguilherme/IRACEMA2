function [R, dR, Jmod] = FastShape3D(GeometryObject,IntegrationPoint, ... 
    global_basis_index, element_local_mapping,element_ranges, element)
qu = IntegrationPoint(1);
qv = IntegrationPoint(2);
qw = IntegrationPoint(3);

pu = GeometryObject.p(1);
pv = GeometryObject.p(2);
pw = GeometryObject.p(3);

U = GeometryObject.knots{1};
V = GeometryObject.knots{2};
W = GeometryObject.knots{3};

support = global_basis_index(element_local_mapping(:,element),:);

u_range = element_ranges(element,:,1);
v_range = element_ranges(element,:,2);
w_range = element_ranges(element,:,3);

u = ((u_range(2)-u_range(1))*qu +(sum(u_range)))/2; % Parent -> Parametric
v = ((v_range(2)-v_range(1))*qv +(sum(v_range)))/2;
w = ((w_range(2)-w_range(1))*qw +(sum(w_range)))/2;

su = FindSpanLinear(length(U)-pu-2,pu,u,U);
sv = FindSpanLinear(length(V)-pv-2,pv,v,V);
sw = FindSpanLinear(length(W)-pw-2,pw,w,W);

P = GeometryObject.points;
ind = sub2ind(size(P),support(:,1),support(:,2),support(:,3));
ActivePoints = P(ind);
ActivePoints = cell2mat(ActivePoints);
Weights = ActivePoints(:,4);
P = ActivePoints(:,1:3);

N = DersBasisFun(su,u,pu,1,U);
M = DersBasisFun(sv,v,pv,1,V);
L = DersBasisFun(sw,w,pw,1,W);

B = kron(kron(L(1,:),M(1,:)),N(1,:));
dBdu = kron(kron(L(1,:),M(1,:)),N(2,:));
dBdv = kron(kron(L(1,:),M(2,:)),N(1,:));
dBdw = kron(kron(L(2,:),M(1,:)),N(1,:));

Q = B*Weights;
dQdu = dBdu*Weights;
dQdv = dBdv*Weights;
dQdw = dBdw*Weights;

R = B'.*Weights/Q;

ratios = Weights/(Q*Q);
dRdu = ratios.*(Q*dBdu' -B'*dQdu);
dRdv = ratios.*(Q*dBdv' -B'*dQdv);
dRdw = ratios.*(Q*dBdw' -B'*dQdw);

x = sum((R.*P));

dxdu = sum(P.*dRdu);
dxdv = sum(P.*dRdv);
dxdw = sum(P.*dRdw);

dXdU = [dxdu', dxdv', dxdw'];
dUdX = inv(dXdU);

dudx = dUdX(:,1)';
dvdx = dUdX(:,2)';
dwdx = dUdX(:,3)';

dR = dRdu*dudx +dRdv*dvdx +dRdw*dwdx;

tmp = element_ranges(element,2,:) - element_ranges(element,1,:);
tmp = squeeze(tmp);
dQdU = 0.5*diag(tmp);
Jacobian = dXdU(:,1)*dQdU(1,:) + dXdU(:,2)*dQdU(2,:) +dXdU(:,3)*dQdU(3,:);
Jmod = det(Jacobian);
end