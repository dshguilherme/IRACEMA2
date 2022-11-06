classdef Assembler
    
    properties
        
        quadrature;
        dimensions;
        domain;
        
    end
     
    methods
        
        function obj = Assembler(quadrature,dimensions,domain)
            
            qrules = ["gauss", "cg"];
            assert(any(strcmp(quadrature,qrules)),"Error: must select a valid quadrature scheme. Valid values are 'gauss' and 'cg'");
            obj.quadrature = quadrature;

            assert(dimensions > 0,"Error: the solutions' dimension must be greater than 0.");
            obj.dimensions = dimensions;

            obj.domain = domain;
        end
        
        function [qpoints, qweights] = quad_rule(obj)
           
            qpoints = cell(obj.domain.rank,1);
            qweights = qpoints;
            
            switch obj.quadrature
                
                case "gauss"
                    for i=1:obj.domain.rank
                        [qp, qw] = gaussian_quadrature(obj.domain.p(i));
                        qpcell{i,1} = qp;
                        qpcell{i,2} = qw;
                    end
                    
                case "cg"
                    for i=1:obj.domain.rank
                       [qp, qw] = cauchy_galerkin_quadrature(obj.domain.p(i));
                        qpcell{i,1} = qp;
                        qpcell{i,2} = qw;
                    end
            end
            
            switch obj.domain.rank
                case 1
                    
                    qpoints = qpcell{1,1};
                    qpoints = qpoints(:);
                    qweights = qpcell{1,2};
                    qweights = qweights(:);
                    
                case 2
                    
                    qu = qpcell{1,1};
                    wu = qpcell{1,2};
                    qv = qpcell{2,1};
                    wv = qpcell{2,2};
                    
                    n_quad = zeros(2,1);
                    n_quad(1) = length(qu);
                    n_quad(2) = length(qv);
                    n_quad = prod(n_quad);
                    
                    qpoints = zeros(n_quad,2);
                    qweights = zeros(n_quad,1);
                    for n =1:n_quad
                       [i,j] = ind2sub([length(qu),length(qv)],n);
                       qpoints(n,:) = [qu(i),qv(j)];
                       qweights(n) = wu(i)*wv(j);
                    end
                    
                case 3
                    
                    qu = qpcell{1,1};
                    wu = qpcell{1,2};
                    qv = qpcell{2,1};
                    wv = qpcell{2,2};
                    qw = qpcell{3,1};
                    ww = qpcell{3,2};
                    
                    n_quad = zeros(3,1);
                    n_quad(1) = length(qu);
                    n_quad(2) = length(qv);
                    n_quad(3) = length(qw);
                    n_quad = prod(n_quad);
                    
                    qpoints = zeros(n_quad,3);
                    qweights = zeros(n_quad,1);
                    
                    for n=1:n_quad
                        [i,j,k] = ind2sub([length(qu),length(qv),length(qw)],n);
                        qpoints(n,:) = [qu(i),qv(j),qw(k)];
                        qweights(n) = wu(i)*wv(j)*ww(k);
                    end
            end
            
        end
        
        function [qpoints, qweights] = boundary_quad_rule(obj)
            qpoints = cell(obj.domain.rank*2,1);
            qweights = qpoints;
            
            for i=1:obj.domain.rank
                [qp, qw] = gaussian_quadrature(2*obj.domain.p(i) +1);
                qpcell{i,1} = qp;
                qpcell{i,2} = qw;
            end
            
            switch obj.domain.rank
                case 1
                    error("For boundaries of 1D geometries, use a Dirichlet Boundary Condition");
                
                case 2
                    qpoints = cell(4,1);
                    qweights = cell(4,1);
                    qu = qpcell{1,1};
                    wu = qpcell{1,2};
                    qv = qpcell{2,1};
                    wv = qpcell{2,2};
                    
                    b1 = ones(size(qu));
                    b2 = ones(size(qu));
                    qpoints{3} = [qu, -b1];                   
                    qweights{3} = wu;
                    qpoints{4} = [qu, b2];
                    qweights{4} = wu;
                    
                    
                    b1 = ones(size(qv));
                    b2 = ones(size(qv));
                    qpoints{1} = [-b1, qv];
                    qweights{1} = wv;
                    qpoints{2} = [b2, qv];
                    qweights{2} = wv;
                
                case 3                   
                    qpoints = cell(6,1);
                    qweights = cell(6,1);
                    qu = qpcell{1,1};
                    wu = qpcell{1,2};
                    qv = qpcell{2,1};
                    wv = qpcell{2,2};
                    qw = qpcell{3,1};
                    ww = qpcell{3,2};
                    
                    nquad = length(qv)*length(qw);
                    tmp1 = zeros(nquad,3);
                    tmp2 = tmp1;
                    tmp3 = zeros(nquad,1);
                    for n=1:n_quad
                        [i,j] = ind2sub([length(qv),length(qw)],n);
                        tmp1(n,:) = [-1 qv(i) qw(j)];
                        tmp2(n,:) = [1 qv(i) qw(j)];
                        tmp3(n) = wv(i)*ww(j);
                    end
                    qpoints{1} = tmp1;
                    qpoints{2} = tmp2;
                    qweights{1} = tmp3;
                    qweights{2} = tmp3;
                    
                    nquad = length(qu)*length(qw);
                    tmp1 = zeros(nquad,3);
                    tmp2 = tmp1;
                    tmp3 = zeros(nquad,1);
                    for n=1:n_quad
                        [i,j] = ind2sub([length(qu),length(qw)],n);
                        tmp1(n,:) = [qu(i) -1 qw(j)];
                        tmp2(n,:) = [qu(i) 1 qw(j)];
                        tmp3(n) = wu(i)*ww(j);
                    end
                    qpoints{3} = tmp1;
                    qpoints{4} = tmp2;
                    qweights{3} = tmp3;
                    qweights{4} = tmp3;
                    
                  nquad = length(qu)*length(qv);
                    tmp1 = zeros(nquad,3);
                    tmp2 = tmp1;
                    tmp3 = zeros(nquad,1);
                    for n=1:n_quad
                        [i,j] = ind2sub([length(qu),length(qv)],n);
                        tmp1(n,:) = [qu(i) qv(j) -1];
                        tmp2(n,:) = [qu(i) qv(j) 1];
                        tmp3(n) = wu(i)*wv(j);
                    end
                    qpoints{5} = tmp1;
                    qpoints{6} = tmp2;
                    qweights{5} = tmp3;
                    qweights{6} = tmp3;     
            end
        end
                    
            
        function id = id_matrix(obj)
            [global_basis_index, element_local_mapping, element_ranges] = ...
                GetConnectivityArrays(obj.domain);
            d = obj.dimensions;
            [id, ~] = BuildGlobalLocalMatrices(element_local_mapping, d);            
        end
        function lm = location_matrix(obj)
            [global_basis_index, element_local_mapping, element_ranges] = ...
                GetConnectivityArrays(obj.domain);
            d = obj.dimensions;
            [~, lm] = BuildGlobalLocalMatrices(element_local_mapping, d);
        end
        
        function F = force_vector(obj, fun)
            d = obj.dimensions;
            [global_basis_index, element_local_mapping, element_ranges] = ...
                GetConnectivityArrays(obj.domain);
            [~, lm] = BuildGlobalLocalMatrices(element_local_mapping, d);
            
            [qp, qw] = obj.quad_rule;
            n_quad = length(qw);
            
            ndof = max(max(element_local_mapping))*d;
            [nel_dof, nel] = size(element_local_mapping);
            
            F = zeros(ndof,1);
           for e=1:nel
                F_e = zeros(d,nel_dof);
                for n=1:n_quad
                    q = qp(n,:);
                    for qq = 1:obj.domain.rank
                        uu = element_ranges(e,:,qq);
                        u(qq) = 0.5*((uu(2)-uu(1))*q(qq) +sum(uu));
                    end
                    x = obj.domain.eval_point(u);
                    [R, ~, J] = FastShape(obj.domain, q, global_basis_index, ...
                        element_local_mapping, element_ranges, e);
                    Jmod = abs(J*qw(n));
                    f = fun(x);
                    for dd=1:d
                        F_e(d,:) = F_e(d,:) +Jmod*f(dd)*(R');
                    end
                end
                idx = lm(:,e)';
                F(idx) = F(idx) +F_e(:);
            end
            F = sparse(F);
        end
        
        function [K, F] = weak_dirichlet(obj,g,boundaries,normal,C,gamma)
            %  WARNING: ONLY WORKING FOR 2D OBJECTS ON XY PLANE
            d = obj.dimensions;
            gbi = obj.domain.global_basis_index;
            [elm e_range] = obj.domain.element_local_mapping;
            id = obj.id_matrix;
            lm = obj.location_matrix;
            
            ndof = max(max(elm))*d;
            [nel_dof, nel] = size(elm);
            K = zeros(ndof);
            F = zeros(ndof,1);
            
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
                    K_e = zeros(nel_dof);
                    F_e = zeros(nel_dof,d);
                    for n=1:nquad
                        q = qq(n,:);
                        [R, dR, J] = FastShape(obj.domain, q, gbi, elm, ...
                                                               e_range, e);
                        Jmod = abs(J*ww(n));
                        while isnan(Jmod) || Jmod > 1e5
                            q(uside) = q(uside)+ ((-1)^bside)*1e-5;
                           
                            [R, dR, J] = FastShape(obj.domain, q, gbi, ...
                                                          elm, e_range, e);
                            Jmod = abs(J*ww(n));
                        end
                        for dd =1:obj.domain.rank
                            uu = e_range(e,:,dd);
                            u(dd) = 0.5*((uu(2)-uu(1))*q(dd) +sum(uu));
                        end
                        x = obj.domain.eval_point(u);
                        f = g(x);
                        ndR = dR*normal';
                        N = kron(R',eye(d));
                        K1 = Jmod*gamma*ndR*N;
                        K2 = Jmod*C*(N'*N);
                        F1 = Jmod*gamma*ndR*f;
                        F2 = Jmod*C*N'*f;
                        for dd=1:d
                            F_e(:,dd) = F_e(:,dd) +F1(:,dd) +F2(:,dd);
                        end
                        K_e = K_e +K1 +K2;
                    end
                    idx = lm(:,e)';
                    K(idx,idx) = K(idx,idx) +K_e;
                    F(idx) = F(idx) +F_e(:);
                end
            end
        end
        
        
        function F = neumann_bc(obj,h,boundaries)
            d = obj.dimensions;
            [qp, qw] = obj.quad_rule;
            gbi = obj.domain.global_basis_index;
            [elm e_range] = obj.domain.element_local_mapping;
            id = obj.id_matrix;
            lm = obj.location_matrix;
            
            ndof = max(max(elm))*d;
            [nel_dof, nel] = size(elm);
            F = zeros(ndof,1);
            
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
                    for n=1:nquad
                        q = qq(n,:);
                        [R, ~, J] = FastShape(obj.domain, q, gbi, elm, ...
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
                        f = h(x);
                        for dd=1:d
                            F_e(:,dd) = F_e(:,dd) +Jmod*f(dd)*R;
                        end
                    end
                    idx = lm(:,e)';
                    F(idx) = F(idx) +F_e(:);
                end
            end
        end       
        
        function [K,F] = robin_bc(obj,r,beta,boundaries)
            d = obj.dimensions;
            [qp, qw] = obj.quad_rule;
            gbi = obj.domain.global_basis_index;
            [elm e_range] = obj.domain.element_local_mapping;
            id = obj.id_matrix;
            lm = obj.location_matrix;
            
            ndof = max(max(elm))*d;
            [nel_dof, nel] = size(elm);
            K = zeros(ndof);
            F = zeros(ndof,1);
            
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
                        [R, ~, J] = FastShape(obj.domain, q, gbi, elm, ...
                                                               e_range, e);
                        Jmod = abs(J*ww(n));
                        while isnan(Jmod) || Jmod > 1e5
                            q(uside) = q(uside)+ ((-1)^bside)*1e-5;
                           
                            [R, ~, J] = FastShape(obj.domain, q, gbi, elm, ...
                                                               e_range, e);
                            Jmod = abs(J*qw(n));
                        end
                        
                        for dd=1:d
                            F_e(:,dd) = F_e(:,dd) +Jmod*r(dd)*R;
                        end
                        K_e = K_e +Jmod*beta*(R*R');
                    end
                    idx = lm(:,e)';
                    K(idx,idx) = K(idx,idx) +K_e; 
                    F(idx) = F(idx) +F_e(:);
                end
            end            
        end
                
        function [d, F, solution] = dirichlet_linear_solve(obj,K,F,g,boundaries)
            boundaries = boundaries(:);
            d = zeros(size(F));
            free_dofs = setdiff(1:length(d),boundaries);
            d(boundaries) = g(:);
            F(free_dofs) = F(free_dofs) - K(free_dofs,boundaries)*d(boundaries);
            d(free_dofs) = K(free_dofs,free_dofs)\F(free_dofs);
            asb = obj;
            solution = Solution(asb, d);
            solution.asb = asb;
        end
           
        function M = L2_projector(obj,nd)
           d = nd;
           [global_basis_index, element_local_mapping, element_ranges] = ...
                GetConnectivityArrays(obj.domain);
           [~, lm] = BuildGlobalLocalMatrices(element_local_mapping, d);
            
           [qp, qw] = obj.quad_rule;
           n_quad = length(qw);
            
           ndof = max(max(element_local_mapping))*d;
           [nel_dof, nel] = size(element_local_mapping);
            
           M = zeros(ndof);
            
           for e=1:nel
               M_e = zeros(nel_dof*d);
               for n=1:n_quad
                   q = qp(n,:);
                   [R, ~, J] = FastShape(obj.domain, q, global_basis_index, ...
                       element_local_mapping, element_ranges, e);
                   Jmod = abs(J*qw(n));
                   N = kron(R',eye(d));
                   M_e = M_e +Jmod*(N'*N);
               end
               idx = lm(:,e)';
               M(idx,idx) = M(idx,idx) +M_e;
           end
           M = sparse(M);         
        end
        
        function F = project_function(obj,fun,nd)
            d = nd;
            [global_basis_index, element_local_mapping, element_ranges] = ...
                GetConnectivityArrays(obj.domain);
           [~, lm] = BuildGlobalLocalMatrices(element_local_mapping, d);
            
           [qp, qw] = obj.quad_rule;
           n_quad = length(qw);
            
           ndof = max(max(element_local_mapping))*d;
           [nel_dof, nel] = size(element_local_mapping);
            
           F = zeros(ndof,1);
           for e=1:nel
               F_e = zeros(d, nel_dof);
               for n=1:n_quad
                   q = qp(n,:);
                   [R, ~, J] = FastShape(obj.domain, q, global_basis_index, ...
                       element_local_mapping, element_ranges, e);
                   Jmod = abs(J*qw(n));
                   for qq = 1:obj.domain.rank
                            uu = element_ranges(e,:,qq);
                            u(qq) = 0.5*((uu(2)-uu(1))*q(qq) +sum(uu));
                   end
                   x = obj.domain.eval_point(u);
                   f = fun(x);
                   for dd=1:d
                            F_e(d,:) = F_e(d,:) +Jmod*f(d)*(R');
                   end
                   idx = lm(:,e)';
                   F(idx) = F(idx) +F_e(:);
               end
           end
        end
        
        function F = project_grad(obj,fun)
            d = 3;
            [global_basis_index, element_local_mapping, element_ranges] = ...
                GetConnectivityArrays(obj.domain);
           [~, lm] = BuildGlobalLocalMatrices(element_local_mapping, d);
            
           [qp, qw] = obj.quad_rule;
           n_quad = length(qw);
            
           ndof = max(max(element_local_mapping))*d;
           [nel_dof, nel] = size(element_local_mapping);
            
           F = zeros(ndof,1);
           for e=1:nel
               F_e = zeros(d, nel_dof);
               for n=1:n_quad
                   q = qp(n,:);
                   [R, dR, J] = FastShape(obj.domain, q, global_basis_index, ...
                       element_local_mapping, element_ranges, e);
                   Jmod = abs(J*qw(n));
                   for qq = 1:obj.domain.rank
                            uu = element_ranges(e,:,qq);
                            u(qq) = 0.5*((uu(2)-uu(1))*q(qq) +sum(uu));
                   end
                   x = obj.domain.eval_point(u);
                   f = fun(x);
                   if size(dR,2) < 3
                       ndim = setdiff([1 2 3], 1:1:size(dR,2));
                       dR = [dR zeros(size(dR,1), numel(ndim))];
                   end
                   F_e = F_e +Jmod*eye(d)*f*dR';
                   idx = lm(:,e)';
                   F(idx) = F(idx) +F_e(:);
               end
           end
        end
        
         function M = derivative_projector(obj,nd,xn)
           d = nd;
           [global_basis_index, element_local_mapping, element_ranges] = ...
                GetConnectivityArrays(obj.domain);
           [~, lm] = BuildGlobalLocalMatrices(element_local_mapping, d);
            
           [qp, qw] = obj.quad_rule;
           n_quad = length(qw);
            
           ndof = max(max(element_local_mapping))*d;
           [nel_dof, nel] = size(element_local_mapping);
            
           M = zeros(ndof);
            
           for e=1:nel
               M_e = zeros(nel_dof*d);
               for n=1:n_quad
                   q = qp(n,:);
                   [R, dR, J] = FastShape(obj.domain, q, global_basis_index, ...
                       element_local_mapping, element_ranges, e);
                   Jmod = abs(J*qw(n));
                   N = kron(R',eye(d));
                   M_e = M_e +Jmod*(N'*dR(:,xn)');
               end
               idx = lm(:,e)';
               M(idx,idx) = M(idx,idx) +M_e;
           end
           M = sparse(M);
        end
        
    end
       
        
end

