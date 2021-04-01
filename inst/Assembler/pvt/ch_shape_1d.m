function [dR, M, nabla_M, nabla_mu, div_nabla_R, J] = ...
            ch_shape_1d(domain, integration_point, global_basis_index, ...
            element_local_mapping, element_ranges, element)                            
 
qu = integration_point(1);

pu = domain.p(1);

U = domain.knots{1};

support = global_basis_index(element_local_mapping(:,element),:);

u_range = element_ranges(element,:,1);
u = ((u_range(2)-u_range(1))*qu +(sum(u_range)))/2; % Parent -> Parametric

su = FindSpanLinear(length(U)-pu-1,pu,u,U);

P = GeometryObject.get_point_cell;
ind = sub2ind(size(P),support(:,1));
ActivePoints = P(ind);
ActivePoints = cell2mat(ActivePoints);
Weights = ActivePoints(:,4);
P = ActivePoints(:,1:3);

N = DersBasisFun(su,u,pu,2,U);

B = N(1,:);
dBdu = N(2,:);
d2Bdu2 = N(3,:);

Q = B*Weights;
dQdu = dBdu*Weights;

R = B'.*Weights/Q;
ratios = Weights/(Q*Q);
dRdu = ratios.*(Q*dBdu' -B'*dQdu);

d2Qdu2 = d2Bdu2*Weights;

d2Rdu2 = ratios.*(Q*d2Bdu2'  -B'*d2Qdu2 -2*dBdu'*dQdu ...
                                                    +(2/Q)*B'*(dQdu*dQdu));
dxdu = P.*dRdu;

dudx = 1./dxdu;

dR = dRdu*dudx;

% This will be ugly. Thanks to curse of dimensionality and my brain
% shutting down on not finding a kronecker product / matrix-vector product

d2Rdx2(:,1) = dudx(1)*(d2Rdu2*dudx(1));         
d2Rdx2(:,2) = dudx(2)*(d2Rdu2*dudx(2));        
d2Rdx2(:,3) = dudx(3)*(d2Rdu2*dudx(3));
        
div_nabla_R = sum(d2Rdx2,2);

tmp = element_ranges(element,2,:) - element_ranges(element,1,:);
J = norm(sum(dxdu))*tmp(1);
M = R.*(1-R);
nabla_M = (1-2*R).*(dR);
nabla_mu = ((0.5/(R -R.^2) -2)').*(dR);

end