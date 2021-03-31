classdef Solution < Geometry
    
    properties
        
        dimension;
                
    end
    
    methods
    
        function obj = Solution(rank, knots, points, p, dimension);
            
            obj@Geometry(rank, knots, points, p);
            obj.dimension = dimension;
        
        end
        
        function eval_solution(obj,paramatric_coordinate_array)

        assert(obj.rank == numel(parametric_coordinate_array), "Error: invalid number of parameters");
        % Eval basis at [u,v,w]
        % Basis * obj.points
        end
    end
end