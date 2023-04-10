classdef CoupledElastoDynamic < Assembler & handle
    properties
        young_modulus;
        poisson_ratio;
        stress_tensor;
        rho; % Density
        alpha; % Rayleigh damping factor
        beta; % Rayleigh damping factor
        lambda; % Lamé parameter
        mu; % Lamé parameter
        D0; % Voigt material properties tensor
        
    end
    
    methods
        
        function obj = CoupledElastoDynamic(young, poisson, density, alpha, beta, ...
                quadrature, dimensions, domain)
            obj@Assembler(quadrature, dimensions, domain);
            n = dimensions;
            
            obj.young_modulus = young;
            obj.poisson_ratio = poisson;
            
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
            obj.D0 = D;
            obj.lambda = lambda;
            obj.mu = mu;
            obj.rho = density;
            obj.alpha = alpha;
            obj.beta = beta;            
        end
        
        %% Assembly Functions
        function [K, M, F] = assembleSystem(obj, f)
            [K, M] = obj.assembleMatrices;
            F = obj.assembleForce(f);
        end
        
        function [K, M] = assembleMatrices(obj)
            d = obj.dimensions;
            [gbi, elm, er, lm, qp, qw, n_quad, ndof, nel_dof, nel] = obj.preLoopParser;
            K = zeros(ndof);
            M = K;
            for e=1:nel
                k = zeros(nel_dof*d);
                m = zeros(nel_dof*d);
                for n=1:n_quad
                    q = qp(n,:);
                    [R, dR, J] = FastShape(obj.domain, q, gbi, elm, er, e);
                    Jmod = abs(J*qw(n));
                    B = elementStiffness(dR, d, nel_dof);
                    N = kron(R', eye(d));
                    D = obj.D0;
                    k = k +Jmod*(B'*D)*B;
                    m = m +Jmod*obj.density*(N'*N);
                    
                end
                idx = lm(:,e)';
                K(idx,idx) = K(idx,idx) +k;
                M(idx,idx) = M(idx,idx) +m;
            end
            K = sparse(K);
            M = sparse(M);
        end
        
        function F = assembleForce(obj, f)
            
        end
        %% Utility Functions
    end
end