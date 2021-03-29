classdef WeakDirichlet < Neumann
    
    properties
        penalty;
    end
    
    methods
        
      function obj = WeakDirichlet(penalty,r)
           
            obj@Neumann(r);
            obj.penalty = penalty;
        
        end        
        
        function F = apply(assembler_object)
        
        end
    end
end