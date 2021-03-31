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
        
        function mu = potencial(obj)
            
            k = 1/(2*obj.theta);
            c = obj.c;
            c = cell2mat(c(:));
            mu = k*log(c./(1-c)) +1 -2*c;
            
        end
        
        function nabla_mu(obj)
        
        end
        
        function M = mobility(obj)
            
            D = obj.diffusivity;
            c = obj.c;
            c = cell2mat(c(:));
            
            M = D*c.*(1-c);
        
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