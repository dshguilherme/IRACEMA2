classdef Poisson < Assembler
    
    properties
       
        alpha;
        
    end
    
    methods
        
        function obj = Poisson(alpha,quadrature,dimensions,domain)
            
            obj@Assembler(quadrature,dimensions,domain);
            
            obj.alpha = alpha;            
                        
        end
        
        function K = build_stiffness(varargin)
                        
        end
        
    end
    
end