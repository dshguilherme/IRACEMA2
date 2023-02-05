classdef CH < Assembler
    properties
        mobility;
        lambda;
    end
    
    
    methods
        
        function obj = CH(quadrature, dimensions, domain)
            assert(dimensions==2,'Cahn-Hilliard assembler only implemented for ndimensions = 2');
            
            obj@Assembler(lambda,mobility,quadrature,dimensions,domain);
            obj.mobility = mobility;
            obj.lambda = lambda;
            
        end
        
        function [c, grad_c, lap_c] = basis_functions(obj)
            d = obj.dimensions;
            [global_basis_index, element_local_mapping, element_ranges] = ...
                GetConnectivityArrays(obj.domain);
                        [~, lm] = BuildGlobalLocalMatrices(element_local_mapping, d);
            
            [qp, qw] = obj.quad_rule;
            n_quad = length(qw);
            
            ndof = max(max(element_local_mapping));
            [nel_dof, nel] = size(element_local_mapping);
            
            K = zeros(ndof);
                   for e=1:nel
                K_e = zeros(nel_dof);
                for n=1:n_quad
                    q = qp(n,:);
                    [R, dR, J] = FastShape(obj.domain, q, global_basis_index, ...
                        element_local_mapping, element_ranges, e);
                    Jmod = abs(J*qw(n));
                    K_e = K_e +Jmod*obj.alpha*(dR*dR');
                end
                idx = lm(:,e)';
                K(idx,idx) = K(idx,idx) +K_e;
            end
            K = sparse(K);
        end
        
    end
    