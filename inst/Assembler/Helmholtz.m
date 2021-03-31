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
        
        function M = build_mass(obj)
            
            d = obj.dimensions;
            [global_basis_index, element_local_mapping, element_ranges] = ...
                GetConnectivityArrays(obj.domain);
            [~, lm] = BuildGlobalLocalMatrices(element_local_mapping, d);
            
            [qp, qw] = obj.quad_rule;
            n_quad = length(qw);
            
            ndof = max(max(element_local_mapping))*d;
            [nel_dof, nel] = size(element_local_mapping);
            
            M = zeros(ndof);
            
            for e=1:nel
                M_e = zeros(nel_dof);
                for n=1:n_quad
                    q = qp(n,:);
                    [R, ~, J] = FastShape(obj.domain, q, global_basis_index, ...
                        element_local_mapping, element_ranges, e);
                    Jmod = abs(J*qw(n));
                    N = kron(R',eye(d));
                    M_e = M_e +Jmod*obj.density*(N'*N);
                end
                idx = lm(:,e)';
                M(idx,idx) = M(idx,idx) +M_e;
            end
            M = sparse(M);
        end
        
        function [K, M] = build_matrices(obj)
            
            d = obj.dimensions;
            [global_basis_index, element_local_mapping, element_ranges] = ...
                GetConnectivityArrays(obj.domain);
            [~, lm] = BuildGlobalLocalMatrices(element_local_mapping, d);
            
            [qp, qw] = obj.quad_rule;
            n_quad = length(qw);
            
            ndof = max(max(element_local_mapping))*d;
            [nel_dof, nel] = size(element_local_mapping);
            
            M = zeros(ndof);
            K = zeros(ndof);
            
            for e=1:nel
                M_e = zeros(nel_dof);
                K_e = zeros(nel_dof);
                for n=1:n_quad
                    q = qp(n,:);
                    [R, dR, J] = FastShape(obj.domain, q, global_basis_index, ...
                        element_local_mapping, element_ranges, e);
                    Jmod = abs(J*qw(n));
                    N = kron(R',eye(d));
                    K_e = K_e +Jmod*obj.alpha*(dR*dR');
                    M_e = M_e +Jmod*obj.density*(N'*N);
                end
                idx = lm(:,e)';
                K(idx,idx) = K(idx,idx) +K_e;
                M(idx,idx) = M(idx,idx) +M_e;
            end
            K = sparse(K);
            M = sparse(M);
        end
        
    end
        
end
    
