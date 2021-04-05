classdef Geometry
   
    properties
        rank;
        knots;
        points;
        p;
        n;
    end
    
    methods
        
        function obj = Geometry(rank, knots, points, p)
            
            assert(ismember(rank,[1, 2, 3]), "Error: invalid tensor rank. Allowed ranks are 1, 2 or 3");
            assert(rank == numel(knots(:)), "Error: you must enter a number of Knots equal to the rank.");
            assert(rank == numel(p), "Error: you must enter a number of polynomial degrees equal to the rank.");
            
            obj.rank = rank;
            obj.knots = knots;
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
                    
                    su = FindSpanLinear(nu,pu,u,U);
                    P = obj.get_point_cell;
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
                    su = FindSpanLinear(nu,pu,u,U);
                    
                    v = parametric_coordinate_array(2);
                    nv = obj.n(2);
                    pv = obj.p(2);
                    V = obj.knots{2};
                    sv = FindSpanLinear(nv,pv,v,V);
                    
                    P = obj.get_point_cell;
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
                    su = FindSpanLinear(nu,pu,u,U);
                    
                    v = parametric_coordinate_array(2);
                    nv = obj.n(2);
                    pv = obj.p(2);
                    V = obj.knots{2};
                    sv = FindSpanLinear(nv,pv,v,V);
                    
                    w = parametric_coordinate_array(3);
                    nw = obj.n(3);
                    pw = obj.p(3);
                    W = obj.knots{3}; 
                    sw = FindSpanLinear(nw,pw,w,W);
                                        
                    P = obj.get_point_cell;
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
        
        function extract_boundaries(obj)
            error('Port from IRACEMA 1.0 underway');
        end
        
        function knot_refine(obj)
            
        end
        
        function degree_elevate(obj)
            
        end
        
        function k_refine(obj)
            degree_elevate()
            knot_refine()
        end
        
        function build_derivative(obj)
            error('In development.');
        end
    end
    
end