classdef Elastodynamic < Elastic
    
    methods
        
        function obj = Elastodynamic(young,poisson,quadrature,dimensions,domain)
            obj@Elastic(young,poisson,quadrature,dimensions,domain);
        end
        
        function M =  build_mass(varargin)
            
        end
        
        function integrate_time(varargin)
            error('Not currently implemented.')
        end
        
    end
    
end
