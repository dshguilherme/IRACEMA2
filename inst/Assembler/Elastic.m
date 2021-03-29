classdef Elastic < Assembler
    
    properties
        
        young_modulus;
        poisson_ratio;
        stress_tensor;
        lambda;
        mu;
        
    end
    
    methods
        
        function obj = Elastic(young,poisson,quadrature,dimensions,domain)
            
            obj@Assembler(quadrature,dimensions,domain);
            
            n = dimensions;
            
            lambda = poisson*young/((1+poisson)*(1-2*poisson));
            mu = young/(2*(1+poisson));
            
            tensor = zeros(size(n));
            tensor(:,:) = 2*mu;
            tensor = tensor+(lambda*eye(n));
            
            obj.stress_tensor = tensor;
            obj.young_modulus = young;
            obj.poisson_ratio = poisson;
            obj.lambda = lambda;
            obj.mu = mu;
            
        end
           
        function K = build_stiffness(varargin)
            
        end
        
    end
    
end