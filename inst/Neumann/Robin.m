classdef Robin < Neumann
    
    properties
        beta;
    end
    
    methods
        
        function obj = Robin(beta,r)
           
            obj@Neumann(r);
            obj.beta = beta;
        
        end
        
        function [K, M] = apply(assembler_object)
        
        end
    end
end