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
        
        function K = build_stiffness(obj)
            
            d = obj.dimensions;
            [global_basis_index, element_local_mapping, element_ranges] = ...
                GetConnectivityArrays(obj.domain);
            [~, lm] = BuildGlobalLocalMatrices(element_local_mapping, d);
            
            [qp, qw] = obj.quad_rule;
            n_quad = length(qw);
            
            ndof = max(max(element_local_mapping))*d;
            [nel_dof, nel] = size(element_local_mapping);
            
            K = zeros(ndof);
            
            for e=1:nel
                K_e = zeros(nel_dof);
                for n=1:n_quad
                    q = qp(n,:);
                    [dR, M, nabla_M, nabla_mu, div_nabla_R, J] = ch_shape(obj.domain, q, global_basis_index, ...
                        element_local_mapping, element_ranges, e);
                    Jmod = abs(J*qw(n));
                    
                    first = dR*((M.*nabla_mu +nabla_M.*div_nabla_R)');
                    second = (div_nabla_R') * (M.*div_nabla_R);
                    K_e = K_e +Jmod*(first +second);
                end
                idx = lm(:,e)';
                K(idx,idx) = K(idx,idx) +K_e;
            end
            K = sparse(K);
        end           

        
        function K = predictor(obj)
        
        end
        
        function K = corrector(obj)
        
        end
        
        function alpha_next = adapt_time(obj)
        
        end
        
        function [d, K] = solve_step()
        end
    end
    
end