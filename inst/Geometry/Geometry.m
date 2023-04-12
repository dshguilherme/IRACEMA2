classdef Geometry < handle
   
    properties
        rank; % int
        knots; % [rank x 1] cell
        points; % [nu x nv x nw] cell
        p; % array
        n; % array
    end
    
    methods
        
        function obj = Geometry(rank, knots, points, p)
            
            assert(ismember(rank,[1, 2, 3]), "Error: invalid tensor rank. Allowed ranks are 1, 2 or 3");
            assert(rank == numel(knots(:)), "Error: you must enter a number of Knots equal to the rank.");
            assert(rank == numel(p), "Error: you must enter a number of polynomial degrees equal to the rank.");
            
            obj.rank = rank;
            for i=1:length(knots)
                num = knots{i} - min(knots{i});
                den = max(knots{i}) - min(knots{i});
                obj.knots{i} = num/den;
            end
            obj.p = p;
            
            n = zeros(rank,1);
            
            for i=1:rank
                n(i) = length(knots{i})-p(i)-1;
            end
            
            obj.n = n;            
            mult = prod(n);
            assert(numel(points(:)) == mult, "Error: invalid number of points. A B-Spline has (n-p) points per parametric direction")
            
            obj.points = points;            
        end
        
        function x = eval_point(obj,parametric_coordinate_array)
            
            assert(obj.rank == numel(parametric_coordinate_array),"Error: invalid number of parameters.")
                        
            switch obj.rank
                
                case 1
                    u = parametric_coordinate_array(1);
                    nu = obj.n(1);
                    pu = obj.p(1);
                    U = obj.knots{1};
                    
                    su = FindSpanLinear(nu-1,pu,u,U);
                    P = obj.points;
                    P = P(:);
                    P = cell2mat(P(su-pu+1:su+1));
                    weights = P(:,4);
                    P = P(:,1:3);
                    B = DersBasisFun(su,u,pu,0,U);
                    Q = B*weights;
                    R = B'.*weights/Q;
                    
                    x = sum((R.*P));
                
                case 2
                             
                    u = parametric_coordinate_array(1);
                    nu = obj.n(1);
                    pu = obj.p(1);
                    U = obj.knots{1}; 
                    su = FindSpanLinear(nu-1,pu,u,U); %Book
                    
                    v = parametric_coordinate_array(2);
                    nv = obj.n(2);
                    pv = obj.p(2);
                    V = obj.knots{2};
                    sv = FindSpanLinear(nv-1,pv,v,V);
                    
                    P = obj.points;
                    P = P(su-pu+1:su+1,sv-pv+1:sv+1);
                    P = cell2mat(P(:));
                    weights = P(:,4);
                    P = P(:,1:3);
                    N = DersBasisFun(su,u,pu,0,U);
                    M = DersBasisFun(sv,v,pv,0,V);
                    B = kron(M,N);
                    Q = B*weights;
                    R = B'.*weights/Q;
                    
                    x = sum((R.*P));
                      
                case 3
                    u = parametric_coordinate_array(1);
                    nu = obj.n(1);
                    pu = obj.p(1);
                    U = obj.knots{1}; 
                    su = FindSpanLinear(nu-1,pu,u,U);
                    
                    v = parametric_coordinate_array(2);
                    nv = obj.n(2);
                    pv = obj.p(2);
                    V = obj.knots{2};
                    sv = FindSpanLinear(nv-1,pv,v,V);
                    
                    w = parametric_coordinate_array(3);
                    nw = obj.n(3);
                    pw = obj.p(3);
                    W = obj.knots{3}; 
                    sw = FindSpanLinear(nw-1,pw,w,W);
                                        
                    P = obj.points;
                    P = P(su-pu+1:su+1,sv-pv+1:sv+1,sw-pw+1:sw+1);
                    P = cell2mat(P(:));
                    weights = P(:,4);
                    P = P(:,1:3);
                    N = DersBasisFun(su,u,pu,0,U);
                    M = DersBasisFun(sv,v,pv,0,V);
                    L = DersBasisFun(sw,w,pw,0,W);
                    B = kron(kron(L,M),N);
                    
                    Q = B*weights;
                    R = B'.*weights/Q;
                    
                    x = sum((R.*P));
           
            end
        end
        
        function x = evalPointFromQuadrature(obj, q, er, e)
                    qu = q(1);
                    U = obj.knots{1};
                    pu = obj.p(1);
                    u_range = er(e,:,1);
                    u = ((u_range(2) - u_range(1))*qu +(sum(u_range)))/2;
                    if obj.rank == 1
                        x = obj.eval_point(u);
                    elseif obj.rank > 1
                        qv = q(2);
                        V = obj.knots{2};
                        pv = obj.p(2);
                        v_range = er(e,:,2);
                        v = ((v_range(2) - v_range(1))*qv +(sum(v_range)))/2;
                        x = obj.eval_point([u v]);
                    elseif obj.rank > 2
                        qw = q(3);
                        W = obj.knots{3};
                        pw = obj.p(3);
                        w_range = er(e,:,3);
                        w = ((w_range(2) - w_range(1))*qw +(sum(w_range)))/2;
                        x = obj.eval_points([u v w]);
                    end
        end

function dx = eval_derivative(obj, parametric_coordinate_array)
            assert(obj.rank == numel(parametric_coordinate_array),"Error: invalid number of parameters.")

          switch obj.rank
              case 1
                    U = obj.knots{1}; 
                    su = FindSpanLinear(nu-1,pu,u,U);

                    P = obj.points;
                    P = P(su-pu+1:su+1);
                    P = cell2mat(P(:));
                    Weights = P(:,4);
                    P = P(:,1:3);

                    N = DersBasisFun(su,u,pu,1,U);

                    B = N(1,:);
                    dBdu = N(2,:);

                    Q = B*Weights;
                    dQdu = dBdu*Weights;

                    R = B'.*Weights/Q;

                    ratios = Weights/(Q*Q);
                    dRdu = ratios.*(Q*dBdu' -B'*dQdu);

                    x = sum((R.*P));

                    dxdu = sum(P.*dRdu);
                    dXdU = dxdu';
                    dUdX = pinv(dXdU);

                    dudx = dUdX(:,1)';

                    dR = dRdu*dudx;
                    dx = dR'*P;                       
                case 2
                    U = obj.knots{1};
                    
                    su = FindSpanLinear(nu-1,pu,u,U);
                    P = obj.points;
                    P = P(:);
                    P = cell2mat(P(su-pu+1:su+1));
                    weights = P(:,4);
                    P = P(:,1:3);
                    B = DersBasisFun(su,u,pu,0,U);
                    Q = B*weights;
                    R = B'.*weights/Q;
                case 2
                             
                    P = obj.points;
                    P = P(su-pu+1:su+1,sv-pv+1:sv+1);
                    P = cell2mat(P(:));
                    Weights = P(:,4);
                    P = P(:,1:3);

                    N = DersBasisFun(su,u,pu,1,U);
                    M = DersBasisFun(sv,v,pv,1,V);

                    B = kron(M(1,:),N(1,:));
                    dBdu = kron(M(1,:),N(2,:));
                    dBdv = kron(M(2,:),N(1,:));

                    Q = B*Weights;
                    dQdu = dBdu*Weights;
                    dQdv = dBdv*Weights;

                    R = B'.*Weights/Q;

                    ratios = Weights/(Q*Q);
                    dRdu = ratios.*(Q*dBdu' -B'*dQdu);
                    dRdv = ratios.*(Q*dBdv' -B'*dQdv);

                    x = sum((R.*P));

                    dxdu = sum(P.*dRdu);
                    dxdv = sum(P.*dRdv);

                    dXdU = [dxdu', dxdv'];
                    dUdX = pinv(dXdU);

                    dR = dRdu*dUdX(1,:) +dRdv*dUdX(2,:);
                    dx = dR'*P;           

                case 3
                    u = parametric_coordinate_array(1);
                    nu = obj.n(1);
                    pu = obj.p(1);
                    U = obj.knots{1}; 
                    su = FindSpanLinear(nu-1,pu,u,U);

                    v = parametric_coordinate_array(2);
                    nv = obj.n(2);
                    pv = obj.p(2);
                    V = obj.knots{2};
                    sv = FindSpanLinear(nv-1,pv,v,V);
                    w = parametric_coordinate_array(3);
                    nw = obj.n(3);
                    pw = obj.p(3);
                    W = obj.knots{3}; 
                    sw = FindSpanLinear(nw-1,pw,w,W);

                    P = obj.points;
                    P = P(su-pu+1:su+1,sv-pv+1:sv+1,sw-pw+1:sw+1);
                    P = cell2mat(P(:));
                    Weights = P(:,4);
                    P = P(:,1:3);

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

                    dR = [dRdu, dRdv, dRdw]*dUdX';
                    dx = dR'*P;
            end
                         
                    P = obj.points;
                    P = P(su-pu+1:su+1,sv-pv+1:sv+1,sw-pw+1:sw+1);
                    P = cell2mat(P(:));
                    weights = P(:,4);
                    P = P(:,1:3);
                    N = DersBasisFun(su,u,pu,0,U);
                    M = DersBasisFun(sv,v,pv,0,V);
                    L = DersBasisFun(sw,w,pw,0,W);
                    B = kron(kron(L,M),N);
                    
                    Q = B*weights;
                    R = B'.*weights/Q;
                    R = R';            
                end
            
             
 function du = parameter_derivative(obj, direction, parametric_coordinate_array)
            assert(obj.rank == numel(parametric_coordinate_array),"Error: invalid number of parameters.")

            switch obj.rank
                case 1
                    u = parametric_coordinate_array(1);
                    nu = obj.n(1);
                    pu = obj.p(1);
                    U = obj.knots{1}; 
                    su = FindSpanLinear(nu-1,pu,u,U);

                    P = obj.points;
                    P = P(su-pu+1:su+1);
                    P = cell2mat(P(:));
                    Weights = P(:,4);
                    P = P(:,1:3);

                    N = DersBasisFun(su,u,pu,1,U);

                    B = N(1,:);
                    dBdu = N(2,:);

                    Q = B*Weights;
                    dQdu = dBdu*Weights;

                    R = B'.*Weights/Q;

                    ratios = Weights/(Q*Q);
                    dRdu = ratios.*(Q*dBdu' -B'*dQdu);
                    du = dRdu'*P;
                    
                case 2
                    u = parametric_coordinate_array(1);
                    nu = obj.n(1);
                    pu = obj.p(1);
                    U = obj.knots{1}; 
                    su = FindSpanLinear(nu-1,pu,u,U);

                    v = parametric_coordinate_array(2);
                    nv = obj.n(2);
                    pv = obj.p(2);
                    V = obj.knots{2};
                    sv = FindSpanLinear(nv-1,pv,v,V);

                    P = obj.points;
                    P = P(su-pu+1:su+1,sv-pv+1:sv+1);
                    P = cell2mat(P(:));
                    Weights = P(:,4);
                    P = P(:,1:3);

                    N = DersBasisFun(su,u,pu,1,U);
                    M = DersBasisFun(sv,v,pv,1,V);

                    B = kron(M(1,:),N(1,:));
                    dBdu = kron(M(1,:),N(2,:));
                    dBdv = kron(M(2,:),N(1,:));

                    Q = B*Weights;
                    dQdu = dBdu*Weights;
                    dQdv = dBdv*Weights;

                    R = B'.*Weights/Q;

                    ratios = Weights/(Q*Q);
                    dRdu = ratios.*(Q*dBdu' -B'*dQdu);
                    dRdv = ratios.*(Q*dBdv' -B'*dQdv);

                    if direction == 1
                        du = dRdu'*P;
                    elseif direction == 2
                        du = dRdv'*P;
                    else
                        error('Unsupported direction')
                    end
                    
                    
                case 3
                    u = parametric_coordinate_array(1);
                    nu = obj.n(1);
                    pu = obj.p(1);
                    U = obj.knots{1}; 
                    su = FindSpanLinear(nu-1,pu,u,U);

                    v = parametric_coordinate_array(2);
                    nv = obj.n(2);
                    pv = obj.p(2);
                    V = obj.knots{2};
                    sv = FindSpanLinear(nv-1,pv,v,V);

                    w = parametric_coordinate_array(3);
                    nw = obj.n(3);
                    pw = obj.p(3);
                    W = obj.knots{3}; 
                    sw = FindSpanLinear(nw-1,pw,w,W);

                    P = obj.points;
                    P = P(su-pu+1:su+1,sv-pv+1:sv+1,sw-pw+1:sw+1);
                    P = cell2mat(P(:));
                    Weights = P(:,4);
                    P = P(:,1:3);

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

                    if direction == 1
                        du = dRdu'*P;
                    elseif direction == 2
                        du = dRdv'*P;
                    elseif direction == 3
                        du = dRdw'*P;
                    else
                        error('Unsupported direction');
                    end
            end
        end

        
        function boundaries = extract_boundaries(obj)
%             Cell array with columns:
%             Geometry | Elements on Boundary | Boundary # | ControlPoint #
%             If boundary # is odd, u/v/w = 0
%             If boundary # is even, u/v/w = 1
%             Boundaries go in u,v,w order, so
%             1 -> u = 0
%             2 -> u = 1
%             3 -> v = 0
%             4 -> v = 1
%             5 -> w = 0
%             6 -> w = 1
            r = obj.rank;
            boundaries = cell(r*2,4);
            elm = obj.element_local_mapping;
            elements = 1:size(elm,2);
            switch r
                case 1
                    boundaries{1,1} = obj.eval_point(0);
                    boundaries{2,1} = obj.eval_point(1);

                case 2
                    n = obj.n;
                    pu = obj.p(1);
                    U = obj.knots{1};

                    pv = obj.p(2);
                    V = obj.knots{2};

                    P = obj.points;
                    
                    P1 = P(1,:);
                    P2 = P(end,:);
                    
                    boundaries{1,1} = Geometry(r-1,{V},P1,[pv]);
                    boundaries{2,1} = Geometry(r-1,{V},P2,[pv]);
                    
                    [s1 s2] = size(P);
                    idx2 = 1:s2;
                    idx1 = ones(size(idx2));
                    
                    b1_points = sub2ind(size(P),idx1,idx2);
                    b2_points = sub2ind(size(P),s1*idx1,idx2);
                    b1_elements = elements(any(ismember(elm,b1_points,"legacy")))';
                    b2_elements = elements(any(ismember(elm,b2_points,"legacy")))';
                    
                    boundaries{1,2} = b1_elements;
                    boundaries{2,2} = b2_elements;
                    boundaries{1,4} = b1_points;
                    boundaries{2,4} = b2_points;
                    
                    P3 = P(:,1);
                    P3 = P3(:);
                    P4 = P(:,end);
                    P4 = P4(:);

                    boundaries{3,1} = Geometry(r-1,{U},P3,[pu]);
                    boundaries{4,1} = Geometry(r-1,{U},P4,[pu]);
                    
                    idx1 = 1:s1;
                    idx2 = ones(size(idx1));

                    b3_points = sub2ind(size(P),idx1,idx2);
                    b4_points = sub2ind(size(P),idx1,s2*idx2);
                    b3_elements = elements(any(ismember(elm,b3_points,"legacy")))';
                    b4_elements = elements(any(ismember(elm,b4_points,"legacy")))';
                    
                    boundaries{3,2} = b3_elements;
                    boundaries{4,2} = b4_elements;
                    boundaries{3,4} = b3_points;
                    boundaries{4,4} = b4_points;
                   clear b1_idx b2_idx                        
                case 3
                    n = obj.n;
                    pu = obj.p(1);
                    U = obj.knots{1};

                    pv = obj.p(2);
                    V = obj.knots{2};

                    pw = obj.p(3);
                    W = obj.knots{3};
                    
                    P = obj.points;
                    
                    [s1 s2 s3] = size(P);

                    P1 = reshape(P(1,:,:),[n(2),n(3)]);
                    P2 = reshape(P(end,:,:),[n(2),n(3)]);
                
                    boundaries{1,1} = Geometry(r-1,{V,W},P1,[pv,pw]);
                    boundaries{2,1} = Geometry(r-1,{V,W},P2,[pv,pw]);
                    
                    idx1 = 1:s1;
                    idx2 = 1:s2;
                    idx3 = 1:s3;
                    
                    [m23 m32] = ndgrid(idx2,idx3);
                    nidx = ones(numel(m23),1);
                    
                    b1_points = sub2ind(size(P),nidx,m23(:),m32(:));
                    b2_points = sub2ind(size(P),s1*nidx,m23(:),m32(:));
                    b1_elements = elements(any(ismember(elm,b1_points,"legacy")))';
                    b2_elements = elements(any(ismember(elm,b2_points,"legacy")))';
                    boundaries{1,2} = b1_elements;
                    boundaries{1,4} = b1_points;
                    boundaries{2,2} = b2_elements;
                    boundaries{2,4} = b2_points;

                    
                    P3 = reshape(P(:,1,:),[n(1),n(3)]);
                    P4 = reshape(P(:,end,:),[n(1),n(3)]);
                                        
                    boundaries{3,1} = Geometry(r-1,{U,W},P3,[pu,pw]);
                    boundaries{4,1} = Geometry(r-1,{U,W},P4,[pu,pw]);

                   [m13 m31] = ndgrid(idx1,idx3);
                    nidx = ones(numel(m13),1);
                     
                    b1_points = sub2ind(size(P),m13(:),nidx,m31(:));
                    b2_points = sub2ind(size(P),m13(:),s2*nidx,m31(:));
                    b1_elements = elements(any(ismember(elm,b1_points,"legacy")))';
                    b2_elements = elements(any(ismember(elm,b2_points,"legacy")))';
                    boundaries{3,2} = b1_elements;
                    boundaries{3,4} = b1_points;
                    boundaries{4,2} = b2_elements;
                    boundaries{4,4} = b2_points;

                    P5 = reshape(P(:,:,1),[n(1),n(2)]);
                    P6 = reshape(P(:,:,end),[n(1),n(2)]);
                    boundaries{5,1} = Geometry(r-1,{U,V},P5,[pu,pv]);
                    boundaries{6,1} = Geometry(r-1,{U,V},P6,[pu,pv]);
                    
                    [m12 m21] = ndgrid(idx1,idx2);
                    nidx = ones(numel(m12),1);
                    
                    b1_points = sub2ind(size(P),m12(:),m21(:),nidx);
                    b2_points = sub2ind(size(P),m12(:),m21(:),s3*nidx);
                    b1_elements = elements(any(ismember(elm,b1_points,"legacy")))';
                    b2_elements = elements(any(ismember(elm,b2_points,"legacy")))';
                    boundaries{5,2} = b1_elements;
                    boundaries{5,4} = b1_points;
                    boundaries{6,2} = b2_elements;
                    boundaries{6,4} = b2_points;
                    clear b1_idx b2_idx
            end
            for i=1:length(boundaries)
                boundaries{i,3} = i;
            end
        end
        
        function obj = knot_refine(obj,knots_to_add,dir)
            switch obj.rank
                case 1
                    for i=1:length(knots_to_add)
                        xi = knots_to_add(i);
                        [obj.n(dir), obj.knots{dir}, obj.points] = ...
                            KnotInsert(obj.n(dir),obj.p(dir), ...
                                       obj.knots{dir},obj.points,xi);
                    end
                    
                case 2
                    for i=1:length(knots_to_add)
                        xi = knots_to_add(i);
                        [obj.n(dir), obj.knots{dir}, obj.points] = ...
                            SurfaceKnotInsert(obj.n(dir),obj.p(dir), ...
                                       obj.knots{dir},obj.points,xi,dir);
                    end
                 
                case 3
                    for i=1:length(knots_to_add)
                        xi = knots_to_add(i);
                        [obj.n(dir), obj.knots{dir}, obj.points] = ...
                            VolumeKnotInsert(obj.n(dir),obj.p(dir), ...
                                       obj.knots{dir},obj.points,xi,dir);
 
                    end
            end
        end
        
        function obj = degree_elevate(obj, t, dir)
            switch obj.rank
                case 1
                    for i=1:t
                        [obj.n(dir), obj.knots{dir}, obj.points] = ... 
                            DegreeElevate(obj.n(dir),obj.p(dir), ...
                                           obj.knots{dir}, obj.points);
                        obj.p(dir) = obj.p(dir)+1;
                    end

                case 2
                    for i=1:t
                        [obj.n(dir), obj.knots{dir}, obj.points] = ...
                            SurfaceDegreeElevate(obj.n(dir),obj.p(dir), ...
                                                 obj.knots{dir}, obj.points, dir);
                         obj.p(dir) = obj.p(dir)+1;
                    end
                   
                case 3
                    for i=1:t
                        [obj.n(dir), obj.knots{dir}, obj.points] = ...
                            VolumeDegreeElevate(obj.n(dir),obj.p(dir), ...
                                                 obj.knots{dir}, obj.points, dir);
                         obj.p(dir) = obj.p(dir)+1;
                    end
            end
        end
        
        function gbi = global_basis_index(obj)
            gbi = zeros(prod(obj.n),obj.rank);
            for i=1:prod(obj.n)
                switch obj.rank
                    case 1
                        gbi(i,1) = ind2sub(obj.n(:)',i);
                    case 2
                        [gbi(i,1) gbi(i,2)] = ind2sub(obj.n(:)',i);
                    case 3
                        [gbi(i,1) gbi(i,2), gbi(i,3)] = ind2sub(obj.n(:)',i); 
                end
            end
        end
        
        function [elm, e_range] = element_local_mapping(obj)
                 switch obj.rank
                    case 1
                        U = obj.knots{1};
                        nu = obj.n(1);
                        pu = obj.p(1);
                        uU = unique(U);
                        elm = zeros(pu+1,numel(uU)-1);
                        e_range = zeros((numel(uU)-1)*(numel(uU)-1),2);
                        for eu=1:numel(uU)-1                               
                            u = 0.5*(uU(eu) +uU(eu+1));
                            sup_u = FindSpanLinear(nu-1,pu,u,U);
                            sup_eu = sup_u-pu+1:sup_u+1;
                            elm(:,eu) = sup_eu;
                            e_range(eu,1) = uU(eu);
                            e_range(eu,2) = uU(eu+1);
                        end
                            
                    case 2
                        U = obj.knots{1};
                        V = obj.knots{2};
                        nu = obj.n(1);
                        nv = obj.n(2);
                        pu = obj.p(1);
                        pv = obj.p(2);
                        uU = unique(U);
                        uV = unique(V);
                        e = 1;
                       elm = zeros(prod(obj.p+1),(numel(uV)-1)*(numel(uU)-1));
                        e_range = zeros((numel(uV)-1)*(numel(uU)-1),2,2);
                        for ev = 1:numel(uV)-1
                            for eu=1:numel(uU)-1
                                u = 0.5*(uU(eu) +uU(eu+1));
                                v = 0.5*(uV(ev) +uV(ev+1));
                                sup_u = FindSpanLinear(nu-1,pu,u,U);
                                sup_v = FindSpanLinear(nv-1,pv,v,V);
                                sup_eu = sup_u-pu+1:sup_u+1;
                                sup_ev = sup_v-pv+1:sup_v+1;
                                [UU, VV] = ndgrid(sup_eu,sup_ev);
                                basis = [UU(:), VV(:)];
                                elm(:,e) = sub2ind(obj.n',basis(:,1),basis(:,2));                            
                                e_range(e,:,1) = [uU(eu) uU(eu+1)];
                                e_range(e,:,2) = [uV(ev) uV(ev+1)];
                                e = e+1;
                            end
                        end
                     case 3
                        U = obj.knots{1};
                        V = obj.knots{2};
                        W = obj.knots{3};
                        nu = obj.n(1);
                        nv = obj.n(2);
                        nw = obj.n(3);
                        pu = obj.p(1);
                        pv = obj.p(2);
                        pw = obj.p(3);
                        uU = unique(U);
                        uV = unique(V);
                        uW = unique(W);
                        elm = zeros(prod(obj.p+1),(numel(uW)-1)*(numel(uV)-1)*(numel(uU)-1));
                        e_range = zeros((numel(uV)-1)*(numel(uU)-1)*(numel(uW)-1),2,3);
                        e = 1;
                        for ew=1:numel(uW)-1
                            for ev=1:numel(uV)-1
                                for eu=1:numel(uU)-1
                                u = 0.5*(uU(eu) +uU(eu+1));
                                v = 0.5*(uV(ev) +uV(ev+1));
                                w = 0.5*(uW(ew) +uW(ew+1));
                                sup_u = FindSpanLinear(nu-1,pu,u,U);
                                sup_v = FindSpanLinear(nv-1,pv,v,V);
                                sup_w = FindSpanLinear(nw-1,pw,w,W);
                                sup_eu = sup_u-pu+1:sup_u+1;
                                sup_ev = sup_v-pv+1:sup_v+1;
                                sup_ew = sup_w-pw+1:sup_w+1;
                                [UU,VV,WW] = ndgrid(sup_eu,sup_ev,sup_ew);
                                basis = [UU(:) VV(:) WW(:)];
                                elm(:,e) = sub2ind(obj.n', basis(:,1),basis(:,2),basis(:,3));
                                e_range(e,:,1) = [uU(eu) uU(eu+1)];
                                e_range(e,:,2) = [uV(ev) uV(ev+1)];
                                e_range(e,:,3) = [uW(ew) uW(ew+1)];
                                e = e+1;
                                end
                            end
                        end
                end 
            end

        
        function obj = uniform_k_refine(obj,Xi,p)
            obj.degree_elevate(p,1);
            obj.knot_refine(Xi,1);
            if obj.rank > 1
                obj.degree_elevate(p,2);
                obj.knot_refine(Xi,2);
            end
            if obj.rank > 2 
               obj.degree_elevate(p,3);
               obj.knot_refine(Xi,3);
            end
            
        end
        
        function h = plot_geo(obj)
            switch obj.rank
                case 1
                    u = linspace(0,1,100);
                    x = zeros(length(u),1);
                    y = x;
                    z = x;
                    for i=1:length(u)
                        point = obj.eval_point(u(i));
                        x(i) = point(1);
                        y(i) = point(2);
                        z(i) = point(3);
                    end
                    h = plot3(x,y,z,'color','black','LineWidth',2);
                case 2
                    u = linspace(0,1,100);
                    v = u;
                    x = zeros(length(u),length(v));
                    y = x;
                    z = x;
                    for i=1:length(u)
                        for j=1:length(v)
                            point = obj.eval_point([u(i) v(j)]);
                            x(i,j) = point(1);
                            y(i,j) = point(2);
                            z(i,j) = point(3);
                        end
                    end
                    h = surf(x,y,z);
                    set(h,'edgecolor','none','facecolor',[0 0 1],'FaceLighting','phong');
                case 3
                   b = obj.extract_boundaries;
                   b = b(:,1);
                   h = cell(6,1);
                   hold on
                   for k=1:6
                       h{k} = b{k}.plot_geo;
                   end
            end
        end
        
        function [] = plot_mesh(obj)
            div = 50;
            switch obj.rank
                case 1
                    u = unique(obj.knots{1});
                    for i=1:numel(u)-1
                        [x,y,z] = obj.parametric_line(i,i+1,div);
                        plot3(x,y,z,'LineWidth',2);
                        hold all;
                    end
                case 2
                    u = unique(obj.knots{1});
                    v = unique(obj.knots{2});
                    for i=1:numel(u)
                        [x,y,z] = obj.parametric_line([i,1],[i,numel(v)],div);
                        plot3(x,y,z,'Color','black','LineWidth',2);
                        hold all;
                    end
                    for i=1:numel(v)
                        [x,y,z] = obj.parametric_line([1,i],[numel(u),i],div);
                        plot3(x,y,z,'Color','black','LineWidth',2);
                        hold all;
                    end
                case 3
                    %extract_boundaries e fazer o mesmo pra rank 2, mas nao
                    %esta funcionando o extract_boundaries!
                    u = unique(obj.knots{1});
                    v = unique(obj.knots{2});
                    w = unique(obj.knots{3});
                    
                    for i=1:numel(u)
                        for j=1:numel(v)
                            for k=1:numel(w)
                                [x,y,z] = obj.parametric_line([1,j,k],[numel(u),j,k],div);
                                plot3(x,y,z,'Color','black','LineWidth',2);
                                hold all;
                                %
                                [x,y,z] = obj.parametric_line([i,1,k],[i,numel(v),k],div);
                                plot3(x,y,z,'Color','black','LineWidth',2);
                                %
                                [x,y,z] = obj.parametric_line([i,j,1],[i,j,numel(w)],div);
                                plot3(x,y,z,'Color','black','LineWidth',2);
                            end
                        end
                    end
            end
        end
        
        function [] = plot_cpoints(obj)
            switch obj.rank
                case 1
                    P = obj.points;
                    P = P(:);
                    P = cell2mat(P(:));
                    plot3(P(:,1,1),P(:,2,1),P(:,3,1),'-','Color','red','LineWidth',1);
                    hold all;
                    plot3(P(:,1,1),P(:,2,1),P(:,3,1),'o','MarkerFaceColor','red');
                case 2
                    P = obj.points;
                    P = P(:);
                    P = cell2mat(P(:));
                    PX = reshape(P(:,1,1),size(obj.points));
                    PY = reshape(P(:,2,1),size(obj.points));
                    PZ = reshape(P(:,3,1),size(obj.points));
                    plot3(PX,PY,PZ,'-','Color','red','LineWidth',2);
                    hold all;
                    plot3(PX',PY',PZ','-','Color','red','LineWidth',2);
                    plot3(P(:,1,1),P(:,2,1),P(:,3,1),'o','MarkerFaceColor','red');
                case 3
                    P = obj.points;
                    P = P(:);
                    P = cell2mat(P(:));
                    plot3(P(:,1,1),P(:,2,1),P(:,3,1),'o','MarkerFaceColor','red');
                    %fix with extract boundaries
                    %for i =1:6, do the same as in rank 2
            end
        end
        
        function [x,y,z] = parametric_line(obj,from_span,to_span, div)
            %for internal use only
            %from_span = [1,2], to_span = [1,4], draw a line from the
            %v_span = 2 to v_span = 4
            %with u fixed at u_span = 1 (1 to 1)
            switch obj.rank
                case 1
                    u = unique(obj.knots{1}); %get non null spans
                    %since it's a curve (or line in parametric space)
                    %length(from_span) and length(to_span)= 1
                    u_ = linspace(u(from_span),u(to_span),div);
                    x = zeros(div,1);
                    y = x;
                    z = x;
                    for i=1:length(x)
                        point = obj.eval_point(u_(i));
                        x(i) = point(1);
                        y(i) = point(2);
                        z(i) = point(3);
                    end
                    
                case 2
                    u = unique(obj.knots{1});
                    v = unique(obj.knots{2});
                    u_ = linspace(u(from_span(1)),u(to_span(1)),div);
                    v_ = linspace(v(from_span(2)),v(to_span(2)),div);
                    x = zeros(div,1);
                    y = x;
                    z = x;
                    for i=1:length(x)
                        point = obj.eval_point([u_(i),v_(i)]);
                        x(i) = point(1);
                        y(i) = point(2);
                        z(i) = point(3);
                    end
                    
                case 3
                    u = unique(obj.knots{1});
                    v = unique(obj.knots{2});
                    w = unique(obj.knots{3});
                    u_ = linspace(u(from_span(1)),u(to_span(1)),div);
                    v_ = linspace(v(from_span(2)),v(to_span(2)),div);
                    w_ = linspace(w(from_span(3)),w(to_span(3)),div);
                    x = zeros(div,1);
                    y = x;
                    z = x;
                    for i=1:length(x)
                        point = obj.eval_point([u_(i),v_(i),w_(i)]);
                        x(i) = point(1);
                        y(i) = point(2);
                        z(i) = point(3);
                    end
            end
        end
        
                       
        
    end
    
end