classdef Elastodynamic < Elastic
    
    properties
        density;
        speed_of_sound;
    end
    
    methods
        
        function obj = Elastodynamic(rho,young,poisson,quadrature,dimensions,domain)
            obj@Elastic(young,poisson,quadrature,dimensions,domain);
            obj.density = rho;
            obj.speed_of_sound = sqrt(obj.young_modulus/rho);
        end
        
        function M =  build_mass(obj)
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
            
           K = zeros(ndof);
           M = zeros(ndof);
           for e=1:nel
                K_e = zeros(nel_dof*d);
                M_e = K_e;
                for n=1:n_quad
                    q = qp(n,:);
                    [R, dR, J] = FastShape(obj.domain, q, global_basis_index, ...
                        element_local_mapping, element_ranges, e);
                    Jmod = abs(J*qw(n));
                    
                    % There should be a better way to write this
                    B = zeros(d,nel_dof*d);
                    for i=1:d
                        B(i,i:d:end) = dR(:,i);
                    end
                    if d == 2
                        B(3,1:d:end) = dR(:,2);
                        B(3,2:d:end) = dR(:,1);
                    elseif d == 3
                        B(4,2:d:end) = dR(:,3);
                        B(4,3:d:end) = dR(:,2);
                        B(5,1:d:end) = dR(:,3);
                        B(5,3:d:end) = dR(:,1);
                        B(6,1:d:end) = dR(:,2);
                        B(6,2:d:end) = dR(:,1);
                    else
                        error('Could not build stress-strain matrix: obj.dimension should be 2 or 3');
                    end
                    C = obj.D;
                    K_e = K_e +Jmod*B'*C*B;
                    N = kron(R',eye(d));
                    M_e = M_e +Jmod*obj.density*(N'*N);
                end
                idx = lm(:,e)';
                K(idx,idx) = K(idx,idx) +K_e;
                M(idx,idx) = M(idx,idx) +M_e;
           end
           K = sparse(K);
           M = sparse(M);            
        end
        
        function [K_c, M_c] = clamp_boundary(obj,K, M)
            boundaries = obj.domain.extract_boundaries;
            clamp_dofs = boundaries(:,2);
            clamp_dofs = unique(cell2mat(clamp_dofs(:)));        
            
            K_c = full(K);
            M_c = full(M);
            
            K_c(clamp_dofs,:) = [];
            K_c(:,clamp_dofs) = [];
            
            M_c(clamp_dofs,:) = [];
            M_c(:,clamp_dofs) = [];
            K_c = sparse(K_c);
            M_c = sparse(M_c);
        end
        
        function [K_c, M_c] = clamp_dofs(obj,K,M,dofs)
            K_c = full(K);
            M_c = full(M);
            
            K_c(dofs,:) = [];
            K_c(:,dofs) = [];
            
            M_c(dofs,:) = [];
            M_c(:,dofs) = [];  
            K_c = sparse(K_c);
            M_c = sparse(M_c);
        end
         function [vecs, omega, solution_cell] = eigensolve(obj,K,M,n_modes,clamp_dofs)
            [K_c, M_c] = obj.clamp_dofs(K,M,clamp_dofs);
            
            [vecs, omega] = eigs(K_c,M_c,n_modes,'sm');
            omega = sqrt(diag(omega));
            
            [s1, s2] = size(vecs);
            dofs = sort(clamp_dofs,'asc');
            ndof = s1+numel(dofs);
            
            padding = zeros(1, s2);
            for i=1:numel(dofs)
                if dofs(i) == 1
                    vecs = [padding; vecs(1:end,:)];
                elseif dofs(i) == ndof
                    vecs = [vecs(1:end,:); padding];
                else
                    vecs = [vecs(1:dofs(i)-1,:); padding; vecs(dofs(i):end,:)];
                end
            end
            solution_cell = cell(n_modes);
            for i=1:n_modes
                solution_cell{i} = Solution(obj, vecs(:,i));
            end
        end          
        function integrate_time(varargin)
            error('Not currently implemented.')
        end
        
    end
    
end
