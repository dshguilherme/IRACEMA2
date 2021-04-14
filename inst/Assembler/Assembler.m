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
        
        function id = id_matrix(obj)
            [global_basis_index, element_local_mapping, element_ranges] = ...
                GetConnectivityArrays(obj.domain);
            d = obj.dimensions;
            [id, ~] = BuildGlobalLocalMatrices(element_local_mapping, d);            
        end
        
        function F = variable_force(obj, fun)
            d = obj.dimensions;
            [global_basis_index, element_local_mapping, element_ranges] = ...
                GetConnectivityArrays(obj.domain);
            [~, lm] = BuildGlobalLocalMatrices(element_local_mapping, d);
            
            [qp, qw] = obj.quad_rule;
            n_quad = length(qw);
            
            ndof = max(max(element_local_mapping))*d;
            [nel_dof, nel] = size(element_local_mapping);
            
            F = zeros(ndof*d,1);
           for e=1:nel
                F_e = zeros(nel_dof,d);
                for n=1:n_quad
                    q = qp(n,:);
                    u = q/2 +0.5;
                    x = obj.domain.eval_point(u);
                    [R, ~, J] = FastShape(obj.domain, q, global_basis_index, ...
                        element_local_mapping, element_ranges, e);
                    Jmod = abs(J*qw(n));
                    for dd=1:d
                        F_e(:,d) = F_e(:,d) +Jmod*fun(x)*(R');
                    end
                end
                idx = lm(:,e)';
                F(idx) = F(idx) +F_e(:);
            end
            F = sparse(F);
        end
        
        function F = constant_force(obj, f)
            d = obj.dimensions;
            assert(length(f) == d,"ERROR: Force should have the same number of dimensions than the solution");
            [global_basis_index, element_local_mapping, element_ranges] = ...
                GetConnectivityArrays(obj.domain);
            [~, lm] = BuildGlobalLocalMatrices(element_local_mapping, d);
            
            [qp, qw] = obj.quad_rule;
            n_quad = length(qw);
            
            ndof = max(max(element_local_mapping))*d;
            [nel_dof, nel] = size(element_local_mapping);
            
            F = zeros(ndof*d,1);
           for e=1:nel
                F_e = zeros(nel_dof,d);
                for n=1:n_quad
                    q = qp(n,:);
                    [R, ~, J] = FastShape(obj.domain, q, global_basis_index, ...
                        element_local_mapping, element_ranges, e);
                    Jmod = abs(J*qw(n));
                    for dd=1:d
                        F_e(:,d) = F_e(:,d) +Jmod*f(d)*(R');
                    end
                end
                idx = lm(:,e)';
                F(idx) = F(idx) +F_e(:);
            end
            F = sparse(F);
        end
        
        function d = L2_projection(obj,F)
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
                    M_e = M_e +Jmod*(R'*R);
                end
                idx = lm(:,e)';
                M(idx,idx) = M(idx,idx) +M_e;
            end
            M = sparse(M);
            d = M\F;
        end
        
        function [d, F, solution] = dirichlet_linear_solve(obj,K,F,g,boundaries)
            d = zeros(size(F));
            free_dofs = setdiff(1:length(d),boundaries);
            d(boundaries) = g;
            F(free_dofs) = F(free_dofs) - K(free_dofs,boundaries)*d(boundaries);
            d(free_dofs) = K(free_dofs,free_dofs)\F(free_dofs);
            solution = Solution(domain, obj.dimensions, d);
        end
       
        
    end
end
