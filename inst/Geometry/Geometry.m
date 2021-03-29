classdef Geometry
   
    properties
        rank;
        knots;
        points;
        p;
    end
    
    methods
        
        function obj = Geometry(rank, knots, points, p)
            
            assert(ismember(rank,[1, 2, 3]), "Error: invalid tensor rank. Allowed ranks are 1, 2 or 3");
            assert(rank == numel(knots(:)), "Error: you must enter a number of Knots equal to the rank.");
            assert(rank == numel(p), "Error: you must enter a number of polynomial degrees equal to the rank.");
            
            obj.rank = rank;
            obj.knots = knots;
            obj.p = p;
            
            mult = 1;
            for i=1:rank
                mult = mult*(length(knots{i})-p(i)-1);
            end
            assert(numel(points(:)) == mult, "Error: invalid number of points. A B-Spline has (n-p) points per parametric direction")
            
            obj.points = points;
            
        end
        
        function eval_point()
        end
        
        function reverse_eval()
        end
        
        function extract_boundaries()
        end
        
        function knot_refine()
        end
        
        function degree_elevate()
        end
        
        function k_refine()
            degree_elevate()
            knot_refine()
        end
        
        function build_derivative()
        end
    end
    
end