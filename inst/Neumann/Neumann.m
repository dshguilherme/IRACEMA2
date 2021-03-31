classdef Neumann
    
    properties
        r;
    end
    
    methods
        
        function obj = Neumann(r)
            obj.r = r;
        end
        
        function F = apply(obj,assembler)
            d = assembler.dimensions;
            [global_basis_index, element_local_mapping, element_ranges] = ...
                GetConnectivityArrays(assembler.domain);
            [~, lm] = BuildGlobalLocalMatrices(element_local_mapping, d);
            
            [qp, qw] = assembler.quad_rule;
            n_quad = length(qw);
            
            ndof = max(max(element_local_mapping))*d;
            [nel_dof, nel] = size(element_local_mapping);
            
            F = zeros(ndof,1);
           for e=1:nel
                F_e = zeros(nel_dof,d);
                for n=1:n_quad
                    q = qp(n,:);
                    [R, ~, J] = FastShape(assembler.domain, q, global_basis_index, ...
                        element_local_mapping, element_ranges, e);
                    Jmod = abs(J*qw(n));
                    f = arrayfun(fun,q);
                    F_e = F_e +Jmod*obj.r*R;
                end
                idx = lm(:,e)';
                F(idx) = F(idx) +F_e(:);
            end
            F = sparse(F);            
        end
        
    end
    
end