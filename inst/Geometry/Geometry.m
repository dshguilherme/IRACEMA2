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
        
        function [x,y,z] = eval_point(parametric_coordinate_array)
            
            assert(obj.rank == numel(parametric_coordinate_array),"Error: invalid number of parameters.")
                        
            switch obj.rank
                
                case 1
                    Pw = [];
                    u = parametric_coordinate_array(1);
                    nu = obj.n(1);
                    pu = obj.p(1);
                    U = obj.knots{1};
                    
                    point = CCurvePoint2(nu,pu,U,Pw,u);
                    x = point.x;
                    y = point.y;
                    z = point.z;
               
                case 2
                    Pw = [];
                    
                    u = parametric_coordinate_array(1);
                    nu = obj.n(1);
                    pu = obj.p(1);
                    U = obj.knots{1}; 
                    
                    v = parametric_coordinate_array(2);
                    nv = obj.n(2);
                    pv = obj.p(2);
                    V = obj.knots{2};
                    
                    point = SurfacePointRAT3(nu,pu,U,nv,pv,V,Pw,u,v);
                    x = point.x;
                    y = point.y;
                    z = point.z;
                
                case 3
                    Pw = [];
                    u = parametric_coordinate_array(1);
                    nu = obj.n(1);
                    pu = obj.p(1);
                    U = obj.knots{1}; 
                    
                    v = parametric_coordinate_array(2);
                    nv = obj.n(2);
                    pv = obj.p(2);
                    V = obj.knots{2};
                    
                    
                    w = parametric_coordinate_array(3);
                    nw = obj.n(3);
                    pw = obj.p(3);
                    W = obj.knots{3};
                    
                    point = VolumePoint(nu,pu,U,nv,pv,V,nw,pw,W,u,v,w);
                    x = point.x;
                    y = point.y;
                    z = point.z;                  
           
            end
        end
        
        function reverse_eval(varargin)
            error('In development');
        end
        
        function extract_boundaries(varargin)
            error('Port from IRACEMA 1.0 underway');
        end
        
        function knot_refine(varargin)
        end
        
        function degree_elevate(varargin)
        end
        
        function k_refine(varargin)
            degree_elevate()
            knot_refine()
        end
        
        function build_derivative(varargin)
            error('In development.');
        end
    end
    
end