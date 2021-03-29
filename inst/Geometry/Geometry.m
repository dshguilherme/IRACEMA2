classdef Geometry
   
    properties
        kind;
        rank;
        knots;
        points;
        p;
    end
    
    methods
        
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