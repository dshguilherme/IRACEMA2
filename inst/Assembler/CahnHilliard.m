classdef CahnHilliard < Assembler
    
    properties
        
        lambda;
        theta;
        rho_0;
        diffusivity;
        L_0;
        alpha;
        c;
        
    end
    
    methods
        
        function mu(obj)
            
        end
        
        function nabla_mu(obj)
        
        end
        
        function mobility(obj)
        
        end
        
        function nabla_mobility(obj)
        
        end
        
        function rho = initial_condition(obj)
        
                        
        end
        
        function K = build_stiffness(obj)
            
        end
        
        function K = predictor(obj)
        
        end
        
        function K = corrector(obj)
        
        end
        
        function alpha_next = adapt_time(obj)
        
        end
    end
    
end