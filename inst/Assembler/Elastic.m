classdef Elastic < Assembler
    
    properties
        
        young_modulus;
        poisson_ratio;
        stress_tensor;
        lambda;
        mu;
        
    end
    
    methods
        
        function obj = Elastic(young,poisson,quadrature,dimensions,domain)
            
            obj@Assembler(quadrature,dimensions,domain);
            
            n = dimensions;
            
            lambda = poisson*young/((1+poisson)*(1-2*poisson));
            mu = young/(2*(1+poisson));
            
            tensor = zeros(2*n -(1-mod(n,2))); % Gambiarra para 2D e 3D
           
            tensor(1:n,1:n) = lambda+2*mu*eye(n);
            tensor(n+1:end,n+1:end) = mu*eye(n -(1-mod(n,2)));
            
            obj.stress_tensor = tensor;
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
                K_e = zeros(nel_dof);
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
                    C = obj.tensor;
                    K_e = K_e +Jmod*B'*C*B;
                end
                idx = lm(:,e)';
                K(idx,idx) = K(idx,idx) +K_e;
           end
           K = sparse(K);
        end
        
    end
    
end