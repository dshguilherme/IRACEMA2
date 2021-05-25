classdef Elastic < Assembler
    
    properties
        
        young_modulus;
        poisson_ratio;
        stress_tensor;
        D;
        lambda;
        mu;
        
    end
    
    methods
        
        function obj = Elastic(young,poisson,quadrature,dimensions,domain)
            
            obj@Assembler(quadrature,dimensions,domain);
            
            n = dimensions;
            
            lambda = poisson*young/((1+poisson)*(1-2*poisson));
            mu = young/(2*(1+poisson));
            
            c = zeros(n,n,n,n);
            D = zeros(2*n -(1-mod(n,2)));
            for i=1:2
                for j=1:2
                    for k=1:2
                        for ell = 1:2
                            c(i,j,k,ell) = ...
                                kronDelta(i,j)*kronDelta(k,ell)*lambda ...
                                +mu*(kronDelta(i,k)*kronDelta(j,ell) ...
                                + kronDelta(i,ell)*kronDelta(j,k));
                            if n == 2
                                I = Voigt2D(i,j);
                                J = Voigt2D(k,ell);
                            elseif n == 3
                                I = Voigt3D(i,j);
                                J = Voigt3D(k,ell);
                            end
                            D(I,J) = c(i,j,k,ell);
                        end
                    end
                end
            end
                   
            obj.stress_tensor = c;
            obj.D = D;
            obj.young_modulus = young;
            obj.poisson_ratio = poisson;
            obj.lambda = lambda;
            obj.mu = mu;
            
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
                K_e = zeros(nel_dof*d);
                for n=1:n_quad
                    q = qp(n,:);
                    [~, dR, J] = FastShape(obj.domain, q, global_basis_index, ...
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
                    D = obj.D;
                    K_e = K_e +Jmod*(B'*D)*B;
                end
                idx = lm(:,e)';
                K(idx,idx) = K(idx,idx) +K_e;
           end
           K = sparse(K);
        end

        function F = variable_neumann_bc(obj,F,h,boundaries)
            d = obj.dimensions;
            asb = Assembler("gauss",d,boundaries{1});
            [qp, qw] = asb.quad_rule;
                 [global_basis_index, element_local_mapping, element_ranges] = ...
                    GetConnectivityArrays(boundaries{1});
                [~, lm] = BuildGlobalLocalMatrices(element_local_mapping,d);
                
                n_quad = length(qw);
                ndof = max(max(element_local_mapping))*d;
                [nel_dof, nel] = size(element_local_mapping);
                Fi = zeros(ndof,1);
                for e=1:nel
                    F_e = zeros(nel_dof,d);
                    for n=1:n_quad
                        q = qp(n,:);
                        [R, ~, J] = FastShape(boundaries{1}, q, global_basis_index, ...
                            element_local_mapping, element_ranges, e);
                        Jmod = abs(J*qw(n));
                        if asb.domain.rank == 1
                            uu = element_ranges(e,:);
                            u = 0.5*((uu(2)-uu(1))*q +sum(uu));
                        else
                            for qq = 1:obj.domain.rank
                                uu = element_ranges(e,:,qq);
                                u(qq) = 0.5*((uu(2)-uu(1))*q(qq) +sum(uu));
                            end
                        end
                        x = asb.domain.eval_point(u);
                        f = h(x);
                        for dd=1:d
                            F_e(:,dd) = F_e(:,dd) + Jmod*f(dd)*R;
                        end
                    end
                    idx = lm(:,e)';
                    Fi(idx) = Fi(idx) +F_e(:);
                end
                id = obj.id_matrix;
                dofs = id(boundaries{2},:);
                F(dofs(:)) = F(dofs(:))+Fi(:);
                F = sparse(F);
        end        
        
        function F = constant_neumann_bc(obj,F,h,boundaries)
            domains = boundaries(:,1);
            d = obj.dimensions;
            assert(length(h) == d,"ERROR: Boundary Condition must have the same size as solution's dimension");

            for i=1:numel(domains)
                asb = Assembler("gauss",d,domains{i});
                [qp, qw] = asb.quad_rule;
                [global_basis_index, element_local_mapping, element_ranges] = ...
                    GetConnectivityArrays(domains{i});
                [~, lm] = BuildGlobalLocalMatrices(element_local_mapping,d);
                
                n_quad = length(qw);
                ndof = max(max(element_local_mapping))*d;
                [nel_dof, nel] = size(element_local_mapping);
                Fi = zeros(d*ndof,1);
                
                for e=1:nel
                    F_e = zeros(nel_dof,d);
                    for n=1:n_quad
                        q = qp(n,:);
                        [R, ~, J] = FastShape(domains{i}, q, global_basis_index, ...
                            element_local_mapping, element_ranges, e);
                        Jmod = abs(J*qw(n));
                        for dd=1:d
                            F_e(:,dd) = F_e(:,dd) +Jmod*h(dd)*(R');
                        end

                    end
                    idx = lm(:,e)';
                    Fi(idx) = Fi(idx) +F_e(:);
                end
                
                id = obj.id_matrix;
                dofs = id(basis);
                dofs = dofs(:);
                F(dofs) = F(dofs)+Fi;
            end
                F = sparse(F);
        end
        
        function [stress, sd] = calculate_stresses(obj, d_solution)
          d = obj.dimensions;
          [global_basis_index, element_local_mapping, element_ranges] = ...
               GetConnectivityArrays(obj.domain);
          [~, lm] = BuildGlobalLocalMatrices(element_local_mapping, d);
          [~, lm2] = BuildGlobalLocalMatrices(element_local_mapping, length(obj.D));
           
           [qp, qw] = obj.quad_rule;
           n_quad = length(qw);
           
           ndof = max(max(element_local_mapping))*d;
           [nel_dof, nel] = size(element_local_mapping);
           
           M = zeros(ndof*length(obj.D)/d);
           F = zeros(ndof*length(obj.D)/d,1);
           stress = zeros(nel,n_quad,length(obj.D));
           for e=1:nel
               Fs_e = zeros(length(obj.D),nel_dof);
               M_e = zeros(nel_dof*length(obj.D));
               for n=1:n_quad
                    q = qp(n,:);
                    [R, dR, J] = FastShape(obj.domain, q, global_basis_index, ...
                        element_local_mapping, element_ranges, e);
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
                    strain = B*d_solution(lm(:,e))*J;
                    stress(e,n,:) = C*strain;
                    for dd=1:length(C)
                        Fs_e(dd,:) = Fs_e(dd,:) +qw(n)*R'*stress(e,n,dd);
                    end
                    N = kron(R',eye(length(C)));
                    M_e = M_e +J*qw(n)*N'*N;
               end
               idx = lm2(:,e)';
               M(idx,idx) = M(idx,idx) +M_e;
               F(idx) = F(idx) +Fs_e(:);
           end
        sd = M\F;
        
    end
    end
end