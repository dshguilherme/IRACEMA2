function [dR, M, nabla_M, nabla_mu, div_nabla_R, J] = ...
            ch_shape_3d(domain, integration_point, global_basis_index, ...
            element_local_mapping, element_ranges, element);
                            
 
qu = integration_point(1);
qv = integration_point(2);
qw = integration_point(3);

pu = domain.p(1);
pv = domain.p(2);
pw = domain.p(3);

U = domain.knots{1};
V = domain.knots{2};
W = domain.knots{3};

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

P = GeometryObject.get_point_cell;
ind = sub2ind(size(P),support(:,1),support(:,2),support(:,3));
ActivePoints = P(ind);
ActivePoints = cell2mat(ActivePoints);
Weights = ActivePoints(:,4);
P = ActivePoints(:,1:3);

N = DersBasisFun(su,u,pu,2,U);
M = DersBasisFun(sv,v,pv,2,V);
L = DersBasisFun(sw,w,pw,2,W);

B = kron(L(1,:),kron(M(1,:),N(1,:)));
dBdu = kron(L(1,:),kron(M(1,:), N(2,:)));
d2Bdu2 = kron(L(1,:),kron(M(1,:), N(3,:)));
d2Bdudv = kron(L(1,:),kron(M(2,:), N(2,:)));
d2Bdudw = kron(L(2,:),kron(M(1,:), N(2,:)));

dBdv = kron(L(1,:),kron(M(2,:), N(1,:)));
d2Bdv2 = kron(L(1,:),kron(M(3,:), N(1,:)));
d2Bdvdw = kron(L(2,:),kron(M(2,:), N(1,:)));

dBdw = kron(L(2,:),kron(M(1,:),N(1,:)));
d2Bdw2 = kron(L(3,:),kron(M(1,:), N(1,:)));

Q = B*Weights;
dQdu = dBdu*Weights;
dQdv = dBdv*Weights;
dQdw = dBdw*Weights;

R = B'.*Weights/Q;
ratios = Weights/(Q*Q);
dRdu = ratios.*(Q*dBdu' -B'*dQdu);
dRdv = ratios.*(Q*dBdv' -B'*dQdv);
dRdw = ratios.*(Q*dBdw' -B'*dQdw);

d2Qdu2 = d2Bdu2*Weights;
d2Qdudv = d2Bdudv*Weights;
d2Qdudw = d2Bdudw*Weights;
d2Qdv2 = d2Bdv2*Weights;
d2Qdvdw = d2Bdvdw*Weights;
d2Qdw2 = d2Bdw2*Weights;

d2Rdu2 = ratios.*(Q*d2Bdu2'  -B'*d2Qdu2 -2*dBdu'*dQdu ...
                                                    +(2/Q)*B'*(dQdu*dQdu));
d2Rdv2 = ratios.*(Q*d2Bdv2'  -B'*d2Qdv2 -2*dBdv'*dQdv ...
                                                    +(2/Q)*B'*(dQdv*dQdv));
d2Rdw2 = ratios.*(Q*d2Bdw2'  -B'*d2Qdw2 -2*dBdw'*dQdw ...
                                                    +(2/Q)*B'*(dQdw*dQdw));
                                                
d2Rdudv = ratios.*(Q*d2Bdudv' -B'*d2Qdudv -dBdu'*dQdv -dBdv'*dQdu ...
                                                    +(2/Q)*B'*(dQdu*dQdv));
                                                
d2Rdudw = ratios.*(Q*d2Bdudw' -B'*d2Qdudw -dBdu'*dQdw -dBdw'*dQdu ...
                                                    +(2/Q)*B'*(dQdu*dQdw));
                                                
d2Rdvdw = ratios.*(Q*d2Bdvdw' -B'*d2Qdvdw -dBdv'*dQdw -dBdw'*dQdv ...
                                                    +(2/Q)*B'*(dQdv*dQdw));
d2Rdvdu = d2Rdudv;
d2Rdwdu= d2Rdudw;
d2Rdwdv = d2Rdvdw;

dxdu = sum(P.*dRdu);
dxdv = sum(P.*dRdv);
dxdw = sum(P.*dRdw);

dXdU = [dxdu', dxdv', dxdw'];
dUdX = inv(dXdU);

dudx = dUdX(:,1)';
dvdx = dUdX(:,2)';
dwdx = dUdX(:,3)';

dR = dRdu*dudx +dRdv*dvdx +dRdw*dwdx;

% This will be ugly. Thanks to curse of dimensionality and my brain
% shutting down on not finding a kronecker product / matrix-vector product

d2Rdx2(:,1) = dudx(1)*(d2Rdu2*dudx(1) +d2Rdudv*dvdx(1) +d2Rdudw*dwdx(1)) + ...
            dvdx(1)*(d2Rdudv*dudx(1) +d2Rdv2*dvdx(1) +d2Rdvdw*dwdx(1)) + ...
            dwdx(1)*(d2Rdudw*dudx(1) +d2Rdvdw*dvdx(1)+d2Rdw2*dwdx(1));
d2Rdx2(:,2) = dudx(2)*(d2Rdu2*dudx(2) +d2Rdudv*dvdx(2) +d2Rdudw*dwdx(2)) + ...
            dvdx(2)*(d2Rdudv*dudx(2) +d2Rdv2*dvdx(2) +d2Rdvdw*dwdx(2)) + ...
            dwdx(2)*(d2Rdudw*dudx(2) +d2Rdvdw*dvdx(2)+d2Rdw2*dwdx(2));
d2Rdx2(:,3) = dudx(3)*(d2Rdu2*dudx(3) +d2Rdudv*dvdx(3) +d2Rdudw*dwdx(3)) + ...
            dvdx(3)*(d2Rdudv*dudx(3) +d2Rdv2*dvdx(3) +d2Rdvdw*dwdx(3)) + ...
            dwdx(3)*(d2Rdudw*dudx(3) +d2Rdvdw*dvdx(3)+d2Rdw2*dwdx(3));
        
%  Won't be needing this, but just in case leave it there for future ref
% d2Rdxdy = dudx(1)*(d2Rdu2*dudx(2) +d2Rdudv*dvdx(2) +d2Rdudw*dwdx(2)) + ...
%            dvdx(1)*(d2Rdudv*dudx(2) +d2Rdv2*dvdx(2) +d2Rdvdw*dwdx(2)) + ...
%            dwdx(1)*(d2Rdudw*dudx(2) +d2Rdvdw*dvdx(2) +d2Rdw2*dwdx(2));
%        
% d2Rdxdz = dudx(1)*(d2Rdu2*dudx(3) +d2Rdudv*dvdx(3) +d2Rdudw*dwdx(3)) + ...
%            dvdx(1)*(d2Rdudv*dudx(3) +d2Rdv2*dvdx(3) +d2Rdvdw*dwdx(3)) + ...
%            dwdx(1)*(d2Rdudw*dudx(3) +d2Rdvdw*dvdx(3) +d2Rdw2*dwdx(3));
%        
% d2Rdydz = dudx(2)*(d2Rdu2*dudx(3) +d2Rdudv*dvdx(3) +d2Rdudw*dwdx(3)) + ...
%            dvdx(2)*(d2Rdudv*dudx(3) +d2Rdv2*dvdx(3) +d2Rdvdw*dwdx(3)) + ...
%            dwdx(2)*(d2Rdudw*dudx(3) +d2Rdvdw*dvdx(3) +d2Rdw2*dwdx(3));
%        
% d2R = [d2Rdx2, d2Rdxdy, d2Rdxdz, d2Rdydz];

div_nabla_R = sum(d2Rdx2,2);

tmp = element_ranges(element,2,:) - element_ranges(element,1,:);
tmp = [squeeze(tmp); 0];
dQdU = eye(3);
dQdU(1,1) = tmp(1);
dQdU(2,2) = tmp(2);
dQdU(3,3) = tmp(3);

Jacobian = dXdU(:,1)*dQdU(1,:) + dXdU(:,2)*dQdU(2,:) +dXdU(:,3)*dQdU(3,:);
J = det(Jacobian);

M = R.*(1-R);
nabla_M = (1-2*R).*(dR);
nabla_mu = ((0.5/(R -R.^2) -2)').*(dR);

end