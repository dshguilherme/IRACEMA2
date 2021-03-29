classdef Helmholtz < Poisson
    
    properties
        density;
        speed_of_sound;
    end
    
    methods
        
        function obj = Helmholtz(rho,alpha,quadrature,dimensions,domain)
            
            obj@Poisson(quadrature,dimensions,domain);
            obj.density = rho;
            obj.speed_of_sound = sqrt(obj.young/rho);
        end
        
        function M = build_mass(varargin)
            
        end
        
    end
    
end
