function [R, dR, Jmod, c, gradc, mu, gradmu] = CHMixedShape(CHAssembler, cmid, mumid, ...
    IntegrationPoint, global_basis_index, element_local_mapping, element_ranges, element)
    
    GeometryObject = CHAssembler.domain;

    qu = IntegrationPoint(1);
    qv = IntegrationPoint(2);

    pu = GeometryObject.p(1);
    pv = GeometryObject.p(2);

    U = GeometryObject.knots{1};
    V = GeometryObject.knots{2};

    support = global_basis_index(element_local_mapping(:,element),:);

    u_range = element_ranges(element,:,1);
    v_range = element_ranges(element,:,2);

    u = ((u_range(2)-u_range(1))*qu +(sum(u_range)))/2; % Parent -> Parametric
    v = ((v_range(2)-v_range(1))*qv +(sum(v_range)))/2;

    su = FindSpanLinear(length(U)-pu-2,pu,u,U);
    sv = FindSpanLinear(length(V)-pv-2,pv,v,V);

    P = GeometryObject.points;
    ind = element_local_mapping(:,element);
    ActivePoints = P(ind);
    ActivePoints = cell2mat(ActivePoints);
    Weights = ActivePoints(:,4);
    P = ActivePoints(:,1:3);

    active_c = cmid(ind);
    active_mu = mumid(length(cmid)/2+ind);

    N = DersBasisFun(su,u,pu,2,U);
    M = DersBasisFun(sv,v,pv,2,V);

    R = kron(M(1,:),N(1,:));
    dRdu = kron(M(1,:),N(2,:));
    dRdv = kron(M(2,:),N(1,:));
    
    c = R*active_c;
    mu = R*active_mu;
    x = R*P;

    dxdu = dRdu*P;
    dxdv = dRdv*P;

    dXdU = [dxdu', dxdv'];
    dUdX = pinv(dXdU);

    dR = dRdu'*dUdX(1,:) +dRdv'*dUdX(2,:);
    dudx = dUdX(1,1);
    dudy = dUdX(1,2);
    dvdx = dUdX(2,1);
    dvdy = dUdX(2,2);
  
    gradc = active_c'*dR;
    gradmu = active_mu'*dR;
    R = R';
    tmp = element_ranges(element,2,:) - element_ranges(element,1,:);
    tmp = squeeze(tmp);
    dQdU = 0.5*diag(tmp);
    Jacobian = dXdU(:,1)*dQdU(1,:) + dXdU(:,2)*dQdU(2,:);
    Jacobian = Jacobian'*Jacobian;
    Jmod = sqrt(det(Jacobian));
end