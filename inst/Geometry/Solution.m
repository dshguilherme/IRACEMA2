classdef Solution < Geometry
    
    properties
        
        dimension;
        d;
                
    end
    
    methods
    
        function obj = Solution(domain, dimension, d)
            
            obj@Geometry(domain.rank, domain.knots, domain.points, domain.p);
            obj.dimension = dimension;
            obj.d = d;
        
        end
        
        function eval_solution(obj,parametric_coordinate_array)

        assert(obj.rank == numel(parametric_coordinate_array), "Error: invalid number of parameters");
        % Eval basis at [u,v,w]
        % Basis * obj.points
        end
    end
end