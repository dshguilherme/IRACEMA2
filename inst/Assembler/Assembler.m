classdef Assembler
    
    properties
        
        quadrature;
        dimensions;
        domain;
        
    end
     
    methods
        
        function obj = Assembler(quadrature,dimensions,domain)
            
            qrules = ["gauss", "cg"];
            assert(isany(quadrature == qrules),"Error: must select a valid quadrature scheme. Valid values are 'gauss' and 'cg'");
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
                        [i,j,k] = ind2sub([length(qu),length(qv),length(qw),n];
                        qpoints(n,:) = [qu(i),qv(j),qw(k)];
                        qweights(n) = wu(i)*wv(j)*ww(k);
                    end
            end
            
        end
              
            
        end
        
        function project(obj, fun)
            
        end
        
        function F = build_force(obj, fun)
            error('Under development.');
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
                F_e = zeros(nel_dof,d);
                for n=1:n_quad
                    q = qp(n,:);
                    [R, ~, J] = FastShape(obj.domain, q, global_basis_index, ...
                        element_local_mapping, element_ranges, e);
                    Jmod = abs(J*qw(n));
                    f = arrayfun(fun,q);
                    F_e = F_e +Jmod*f*(R');
                end
                idx = lm(:,e)';
                F(idx) = F(idx) +F_e(:);
            end
            F = sparse(F);
        end
        
        function apply_dirichlet()
            
        end
        
        function solve()
            
        end
        
    end
    
end