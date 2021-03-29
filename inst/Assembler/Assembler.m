classdef Assembler
    
    properties
        
        quadrature;
        dimensions;
        domain;
        
    end
    
    methods
        
        function obj = Assembler(quadrature,dimensions,domain)
            qrules = ["gauss", "cg"];
            assert(isany(quadrature == qrules),"Error: must select a valid quadrature scheme. Valid values are 'gauss' and 'cg'");
            obj.quadrature = quadrature;

            assert(dimensions > 0,"Error: the solutions' dimension must be greater than 0.");
            obj.dimensions = dimensions;

            obj.domain = domain;
        end
        
        function build_matrices()
        
        end
        
        function project()
            
        end
        
        function build_force()
        
        end
        
        function apply_dirichlet()
            
        end
        
        function solve()
            
        end
        
    end
    
end