classdef Elastodynamic < Elastic
    
    properties
        density;
        speed_of_sound;
    end
    
    methods
        
        function obj = Elastodynamic(rho,young,poisson,quadrature,dimensions,domain)
            obj@Elastic(young,poisson,quadrature,dimensions,domain);
            obj.density = rho;
            obj.speed_of_sound = sqrt(obj.young/rho);
        end
        
        function M =  build_mass(varargin)
            
        end
        
        function integrate_time(varargin)
            error('Not currently implemented.')
        end
        
    end
    
end
