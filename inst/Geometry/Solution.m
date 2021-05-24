classdef Solution < Geometry
    
    properties
        
        id;
        d;
        asb; % Assembler object
        domain; % Geometry obj
        cpoints; % Solution's dofs
    end
    
    methods
    
        function obj = Solution(asb, d)
            
            obj@Geometry(asb.domain.rank, asb.domain.knots, ...
                         asb.domain.points, asb.domain.p);
            obj.domain = asb.domain;
            obj.id = asb.id_matrix;
            obj.d = d;
            
            [s2, ~] = size(obj.id);
            cpoints = cell(1,s2);
            for i=1:s2
                idx = obj.id(i,:);
                cpoints{i} = [d(idx) 1];
            end
            cpoints = reshape(cpoints,size(obj.points));
            obj.cpoints = cpoints;
        end
        
        function [x, d] = eval_solution(obj,parametric_coordinate_array)
            
            x = obj.domain.eval_point(parametric_coordinate_array);
            
            switch obj.rank
                
                case 1
                    u = parametric_coordinate_array(1);
                    nu = obj.n(1);
                    pu = obj.p(1);
                    U = obj.knots{1};
                    
                    su = FindSpanLinear(nu-1,pu,u,U);
                    P = obj.cpoints;
                    P = P(:);
                    P = cell2mat(P(su-pu+1:su+1));
                    weights = P(:,end);
                    P = P(:,1:end-1);
                    B = DersBasisFun(su,u,pu,0,U);
                    Q = B*weights;
                    R = B'.*weights/Q;
                    
                    d = sum((R.*P));  
                
                case 2
                    u = parametric_coordinate_array(1);
                    nu = obj.n(1);
                    pu = obj.p(1);
                    U = obj.knots{1}; 
                    su = FindSpanLinear(nu-1,pu,u,U); %Book
                    
                    v = parametric_coordinate_array(2);
                    nv = obj.n(2);
                    pv = obj.p(2);
                    V = obj.knots{2};
                    sv = FindSpanLinear(nv-1,pv,v,V);
                    
                    P = obj.cpoints;
                    P = P(su-pu+1:su+1,sv-pv+1:sv+1);
                    P = cell2mat(P(:));
                    weights = P(:,end);
                    P = P(:,1:end-1);
                    N = DersBasisFun(su,u,pu,0,U);
                    M = DersBasisFun(sv,v,pv,0,V);
                    B = kron(M,N);
                    Q = B*weights;
                    R = B'.*weights/Q;
                    
                    d = sum((R.*P));               
                
                case 3
                    u = parametric_coordinate_array(1);
                    nu = obj.n(1);
                    pu = obj.p(1);
                    U = obj.knots{1}; 
                    su = FindSpanLinear(nu-1,pu,u,U);
                    
                    v = parametric_coordinate_array(2);
                    nv = obj.n(2);
                    pv = obj.p(2);
                    V = obj.knots{2};
                    sv = FindSpanLinear(nv-1,pv,v,V);
                    
                    w = parametric_coordinate_array(3);
                    nw = obj.n(3);
                    pw = obj.p(3);
                    W = obj.knots{3}; 
                    sw = FindSpanLinear(nw-1,pw,w,W);
                                        
                    P = obj.cpoints;
                    P = P(su-pu+1:su+1,sv-pv+1:sv+1,sw-pw+1:sw+1);
                    P = cell2mat(P(:));
                    weights = P(:,end);
                    P = P(:,1:end-1);
                    N = DersBasisFun(su,u,pu,0,U);
                    M = DersBasisFun(sv,v,pv,0,V);
                    L = DersBasisFun(sw,w,pw,0,W);
                    B = kron(L,kron(M,N));
                    
                    Q = B*weights;
                    R = B'.*weights/Q;
                    
                    d = sum((R.*P));    
            end
        end
        
        function solution_cell = extract_solution_boundaries(obj)
            assert(obj.domain.rank == 3, "This function is only usable for rank 3 solutions");
            b = obj.domain.extract_boundaries;
            dof = b(:,2);
            b = b(:,1);
            d = size(obj.id,1);
            solution_cell = cell(numel(b),1);
            for i=1:numel(b)
                asb = Assembler("gauss",d,b{i});
                dofs = obj.id(:,dof{i});
                s = obj.d(dofs);
                solution_cell{i} = Solution(asb,s);
            end
        end
        
        function h = plot_solution(obj,dim)
           switch obj.rank
                case 1
                    u = linspace(0,1,100);
                    x = zeros(length(u),1);
                    y = x;
                    z = x;
                    for i=1:length(u)
                        point = obj.eval_solution(u(i));
                        x(i) = point(1);
                        y(i) = point(2);
                        z(i) = point(3);
                    end
                    h = plot3(x,y,z,'color','black','LineWidth',2);
                 case 2
                    u = linspace(0,1,100);
                    v = u;
                    x = zeros(length(u),length(v));
                    y = x;
                    z = x;
                    c = x;
                    for i=1:length(u)
                        for j=1:length(v)
                            [point, d] = obj.eval_solution([u(i) v(j)]);
                            x(i,j) = point(1);
                            y(i,j) = point(2);
                            z(i,j) = point(3);
                            c(i,j) = d(dim);
                        end
                    end
                    h = surf(x,y,z,c);
                    colorbar;
                    set(h,'edgecolor','none','FaceLighting','phong');
                
               case 3
                   b_sols = obj.extract_solution_boundaries;
                   hold on
                   for k=1:6
                      h{k} = b_sols{k}.plot_solution(dim);
                   end

           end
        end
        
        function [error_norm, energy_norm] = l2_error_norm(obj,fun,dim)
            asb = obj.asb;
            d = asb.dimensions;
         [global_basis_index, element_local_mapping, element_ranges] = ...
                GetConnectivityArrays(asb.domain);
         [~, lm] = BuildGlobalLocalMatrices(element_local_mapping, d);
         [qp, qw] = asb.quad_rule;
         n_quad = length(qw);
         [~, nel] = size(element_local_mapping);
         error_norm = 0;
         energy_norm = 0;
         for e=1:nel
             for n=1:n_quad
                 q = qp(n,:);
                 csi = q/2 +0.5;
                 [x, uh] = obj.eval_solution(csi);
                 error_norm = qw(n)*(fun(x)-uh(dim))^2;
                 energy_norm = qw(n)*(uh(dim)*uh(dim));
             end
         end
         energy_norm = sqrt(energy_norm);
         error_norm = sqrt(error_norm)/energy_norm;
        end
        
    end
end