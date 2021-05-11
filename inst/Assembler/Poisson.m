classdef Poisson < Assembler
    
    properties
       
        alpha;
        
    end
    
    methods
        
        function obj = Poisson(alpha,quadrature,dimensions,domain)
            
            obj@Assembler(quadrature,dimensions,domain);
            
            obj.alpha = alpha;            
                        
        end
        
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
                    [~, dR, J] = FastShape(obj.domain, q, global_basis_index, ...
                        element_local_mapping, element_ranges, e);
                    Jmod = abs(J*qw(n));
                    K_e = K_e +Jmod*obj.alpha*(dR*dR');
                end
                idx = lm(:,e)';
                K(idx,idx) = K(idx,idx) +K_e;
            end
            K = sparse(K);
        end
        
        function F = neumann_bc(F,obj,h,boundaries)
            domains = boundaries(:,1);
            d = obj.dimensions;
            for i=1:numel(domains)
                asb = Assembler("gauss",d,domains{i});
                [qp, qw] = asb.quad_rule;
                [global_basis_index, element_local_mapping, element_ranges] = ...
                    GetConnectivityArrays(domains{i});
                [~, lm] = BuildGlobalLocalMatrices(element_local_mapping,d);
                
                n_quad = length(qw);
                ndof = max(max(element_local_mapping))*d;
                [nel_dof, nel] = size(element_local_mapping);
                Fi = zeros(ndof,1);
                
                for e=1:nel
                    F_e = zeros(nel_dof,1);
                    for n=1:n_quad
                        q = qp(n,:);
                        [R, ~, J] = FastShape(domains{i}, q, global_basis_index, ...
                            element_local_mapping, element_ranges, e);
                        Jmod = abs(J*qw(n));
                        F_e = F_e +Jmod*h*(R');
                    end
                    idx = lm(:,e)';
                    Fi(idx) = Fi(idx) +F_e(:);
                end
                
                basis = boundaries{i,2};
                basis = basis(:);
                id = obj.id_matrix;
                dofs = id(basis);
                dofs = dofs(:);
                F(dofs) = F(dofs)+Fi;
            end
                F = sparse(F);
        end
            
        function [K,F] = robin_bc(obj,K,F,r,beta,boundaries)
            domains = boundaries(:,1);
            d = obj.dimensions;
            for i=1:numel(domains)
                asb = Assembler("gauss",d,domains{i});
                [qp, qw] = asb.quad_rule;
                [global_basis_index, element_local_mapping, element_ranges] = ...
                    GetConnectivityArrays(domains{i});
                [~, lm] = BuildGlobalLocalMatrices(element_local_mapping,d);

                n_quad = length(qw);
                ndof = max(max(element_local_mapping))*d;
                [nel_dof, nel] = size(element_local_mapping);
                Fi = zeros(ndof*d,1);
                Ki = zeros(ndof);

                for e=1:nel
                    K_e = zeros(nel_dof);
                    F_e = zeros(nel_dof,d);
                    for n=1:n_quad
                        q = qp(n,:);
                        [R, ~, J] = FastShape(domains{i}, q, global_basis_index, ...
                            element_local_mapping, element_ranges, e);
                        Jmod = abs(J*qw(n));
                        for i=1:d
                            F_e(:,d) = F_e(:,d) +Jmod*r(d)*(R);
                        end
                        N = kron(R', eye(d));
                        K_e = K_e +Jmod*beta*(N'*N);
                    end
                    idx = lm(:,e)';
                    Ki(idx,idx) = Ki(idx,idx) + K_e;
                    Fi(idx) = Fi(idx) +F_e(:);
                end

                basis = boundaries{i,2};
                basis = basis(:);
                id = obj.id_matrix;
                dofs = id(basis);
                dofs = dofs(:);                
                F(dofs) = F(dofs)+Fi;
                K(dofs,dofs) = K(dofs,dofs) +Ki;
            end
                F = sparse(F);
                K = sparse(K);
                
            end
        end
        
end