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
                    B = kron(L,kron(M,N));
                    
                    Q = B*weights;
                    R = B'.*weights/Q;
                    
                    x = sum((R.*P));
           
            end
        end
        
        function reverse_eval(obj, physical_coordinate_array)
            error('In development');
%           Step 1. Is point(x,y,z) inside Convex Hull?
%           Step 2. KnotRefine -> Bezier patches
%           Step 3. Find the neareast candidate points
%           Step 4. Find the closest point to the curve within tolerance
        end
        
        function boundaries = extract_boundaries(obj)
            r = obj.rank;
            boundaries = cell(r*2,2);
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

                    P3 = P(:,1);
                    P3 = P3(:);
                    P4 = P(:,end);
                    P4 = P4(:);

                    boundaries{3,1} = Geometry(r-1,{U},P3,[pu]);
                    boundaries{4,1} = Geometry(r-1,{U},P4,[pu]);

                case 3
                    n = obj.n;
                    pu = obj.p(1);
                    U = obj.knots{1};

                    pv = obj.p(2);
                    V = obj.knots{2};

                    pw = obj.p(3);
                    W = obj.knots{3};
                    
                    P = obj.points;
         
                    P1 = reshape(P(1,:,:),[n(2),n(3)]);
                    P2 = reshape(P(end,:,:),[n(2),n(3)]);
                    P3 = reshape(P(:,1,:),[n(1),n(3)]);
                    P4 = reshape(P(:,end,:),[n(1),n(3)]);
                    P5 = reshape(P(:,:,1),[n(1),n(2)]);
                    P6 = reshape(P(:,:,end),[n(1),n(2)]);

                    boundaries{1,1} = Geometry(r-1,{V,W},P1,[pv,pw]);
                    boundaries{2,1} = Geometry(r-1,{V,W},P2,[pv,pw]);
                    boundaries{3,1} = Geometry(r-1,{U,W},P3,[pu,pw]);
                    boundaries{4,1} = Geometry(r-1,{U,W},P4,[pu,pw]);
                    boundaries{5,1} = Geometry(r-1,{U,V},P5,[pu,pv]);
                    boundaries{6,1} = Geometry(r-1,{U,V},P6,[pu,pv]);
            end
            tmp = GetBoundaries(obj);
            boundaries(:,2) = tmp;
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
        
        function obj = uniform_k_refine(obj,Xi,p)
            obj.degree_elevate(p,1);
            obj.degree_elevate(p,2);
            obj.degree_elevate(p,3);
            obj.knot_refine(Xi,1);
            obj.knot_refine(Xi,2);
            obj.knot_refine(Xi,3);
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
                    u = linspace(0,1,100);
                    v = u;
                    w = u;
                    h = cell(6,1);
                    x = zeros(length(u),length(v));
                    y = x;
                    z = x;
                    %%% u,v, w = 0
                    for i=1:length(u)
                        for j=1:length(v)
                            point = obj.eval_point([u(i) v(j) 0]);
                            x(i,j) = point(1);
                            y(i,j) = point(2);
                            z(i,j) = point(3);
                        end
                    end
                    h{1} = surf(x,y,z);
                    hold all;
                    set(h{1},'edgecolor','none','facecolor',[0 0 1],'FaceLighting','phong');
                    %%% u,v, w = 1
                    for i=1:length(u)
                        for j=1:length(v)
                            point = obj.eval_point([u(i) v(j) 1]);
                            x(i,j) = point(1);
                            y(i,j) = point(2);
                            z(i,j) = point(3);
                        end
                    end
                    h{2} = surf(x,y,z);
                    set(h{2},'edgecolor','none','facecolor',[0 0 1],'FaceLighting','phong');
                    %%% u,w, v = 0
                    for i=1:length(u)
                        for j=1:length(w)
                            point = obj.eval_point([u(i) 0 w(j)]);
                            x(i,j) = point(1);
                            y(i,j) = point(2);
                            z(i,j) = point(3);
                        end
                    end
                    h{3} = surf(x,y,z);
                    set(h{3},'edgecolor','none','facecolor',[0 0 1],'FaceLighting','phong');
                    %%% u,w, v = 1
                    for i=1:length(u)
                        for j=1:length(w)
                            point = obj.eval_point([u(i) 1 w(j)]);
                            x(i,j) = point(1);
                            y(i,j) = point(2);
                            z(i,j) = point(3);
                        end
                    end
                    h{4} = surf(x,y,z);
                    set(h{4},'edgecolor','none','facecolor',[0 0 1],'FaceLighting','phong');
                    %%% w,v, u = 0
                    for i=1:length(w)
                        for j=1:length(v)
                            point = obj.eval_point([0 v(j) w(i)]);
                            x(i,j) = point(1);
                            y(i,j) = point(2);
                            z(i,j) = point(3);
                        end
                    end
                    %%% w,v, u = 1
                    h{5} = surf(x,y,z);
                    set(h{5},'edgecolor','none','facecolor',[0 0 1],'FaceLighting','phong');
                    for i=1:length(w)
                        for j=1:length(v)
                            point = obj.eval_point([1 v(j) w(i)]);
                            x(i,j) = point(1);
                            y(i,j) = point(2);
                            z(i,j) = point(3);
                        end
                    end
                    h{6} = surf(x,y,z);
                    set(h{6},'edgecolor','none','facecolor',[0 0 1],'FaceLighting','phong');
                    light;
            end
        end
        
    end
    
end