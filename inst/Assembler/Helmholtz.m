classdef Helmholtz < Poisson
    
    properties
        density;
        speed_of_sound;
    end
    
    methods
        
        function obj = Helmholtz(rho,alpha,quadrature,dimensions,domain)
            
            obj@Poisson(alpha,quadrature,dimensions,domain);
            obj.density = rho;
            obj.speed_of_sound = sqrt(obj.alpha/rho);
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
            vecs = vecs./(max(abs(vecs)));
            for i=1:size(vecs,2)
                if sum(sign(vecs(:,i))) < 0
                    vecs(:,i) = -vecs(:,i);
                end
            end
            solution_cell = cell(n_modes,1);
            for i=1:n_modes
                solution_cell{i} = Solution(obj, vecs(:,i));
            end
        end
        
        function [d, F, solution] = helmholtz_linear_solve(obj,K,M,F,w,g,boundaries)
            boundaries = boundaries(:);
            d = zeros(size(F));
            free_dofs = setdiff(1:length(d),boundaries);
            d(boundaries) = g(:);
            F(free_dofs) = F(free_dofs) - K(free_dofs,boundaries)*d(boundaries);
            d(free_dofs) = (K(free_dofs,free_dofs)  ...
                               -(w^2)*M(free_dofs,free_dofs))\F(free_dofs);
            asb = obj;
            solution = Solution(asb, d);
            solution.asb = asb;
        end
        
        function [K,F] = lui_bc(obj,f,g,normal,boundaries)
            d = obj.dimensions;
            [qp, qw] = obj.quad_rule;
            gbi = obj.domain.global_basis_index;
            [elm e_range] = obj.domain.element_local_mapping;
            id = obj.id_matrix;
            lm = obj.location_matrix;
            
            ndof = max(max(elm))*d;
            [nel_dof, nel] = size(elm);
            F = zeros(ndof,1);
            K = zeros(ndof);
                        
            [nb, ~, ~] = size(boundaries);
            [qp qw] = obj.boundary_quad_rule;           
            for b=1:nb
                elements = boundaries{b,2};
                bside = boundaries{b,3};
                uside = ceil(bside/2);
                qq = qp{bside};
                ww = qw{bside};
                nquad = length(ww);
                for ee=1:numel(elements)
                    e = elements(ee);
                    F_e = zeros(nel_dof,d);
                    K_e = zeros(nel_dof);
                    for n=1:nquad
                        q = qq(n,:);
                        [R, dR, J] = FastShape(obj.domain, q, gbi, elm, ...
                                                               e_range, e);
                        Jmod = abs(J*ww(n));
                        while isnan(Jmod) || Jmod > 1e5
                            q(uside) = q(uside)+ ((-1)^bside)*1e-5;
                           
                            [R, ~, J] = FastShape(obj.domain, q, gbi, elm, ...
                                                               e_range, e);
                            Jmod = abs(J*ww(n));
                        end
                        
                        for dd =1:obj.domain.rank
                            uu = e_range(e,:,dd);
                            u(dd) = 0.5*((uu(2)-uu(1))*q(dd) +sum(uu));
                        end
                        x = obj.domain.eval_point(u);
                        tmp1 = g(x);
                        tmp2 = f(x);
                        r = dot(tmp2,normal) +tmp1;
                        beta = 1;
                        for dd=1:d
                            F_e(:,dd) = F_e(:,dd) +Jmod*1e12*r(dd)*R;
                        end
                        N = kron(R',eye(d));
                        K_e = K_e -Jmod*1e12*N'*N;
                    end
                    idx = lm(:,e)';
                    F(idx) = F(idx) +F_e(:);
                    K(idx,idx) = K(idx,idx) +K_e;
                end
            end          
        end
        
    end
        
end
    
