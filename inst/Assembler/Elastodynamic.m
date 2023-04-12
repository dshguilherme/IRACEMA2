classdef Elastodynamic < Assembler & handle
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
        
        function obj = Elastodynamic(young, poisson, density, alpha, beta, ...
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
                    m = m +Jmod*obj.rho*(N'*N);
                    
                end
                idx = lm(:,e)';
                K(idx,idx) = K(idx,idx) +k;
                M(idx,idx) = M(idx,idx) +m;
            end
            K = sparse(K);
            M = sparse(M);
        end
        
        function F = assembleForce(obj, f)
            d = obj.dimensions;
            [gbi, elm, er, lm, qp, qw, n_quad, ndof, nel_dof, nel] = obj.preLoopParser;
            F = zeros(ndof,1);
            for e=1:nel
                f_e = zeros(nel_dof,d);
                for n=1:n_quad
                    q = qp(n,:);
                    [R, dR, J] = FastShape(obj.domain, q, gbi, elm, er, e);
                    Jmod = abs(J*qw(n));
                    x = obj.domain.evalPointFromQuadrature(q, er, e);
                    ff = f(x);
                    for i=1:d
                        f_e(:,i) = f_e(:,i) +Jmod*R*ff(i);
                    end
                end
                idx = lm(:,e)';
                F(idx) = F(idx) + f_e(:);
            end
            F = sparse(F);
        end
        %% Boundary Conditions
        function F = applyTraction(obj, t, side)
            % Apply Neumann boundary condition on one side of the domain
            d = obj.dimensions;
            [gbi, elm, er, lm, ~, ~, ~, ndof, nel_dof, nel] = obj.preLoopParser;
            F = zeros(ndof,1);
            b = obj.domain.extract_boundaries;
            elements = b{side,2};
            [qp, qw] = obj.boundary_quad_rule(side);
            n_quad = length(qw);
            for i = 1:numel(elements)
                e = elements(i);
                f_e = zeros(nel_dof, d);
                for n=1:n_quad
                    q = qp(n,:);
                    [R, dR, J] = FastShape(obj.domain, q, gbi, elm, er, e);
                    Jmod = abs(J*qw(n));
                    for j=1:d
                        f_e(:,j) = f_e(:,j) +Jmod*R*t(j);
                    end
                end
                idx = lm(:,e)';
                F(idx) = F(idx) +f_e(:);
            end
            F = sparse(F);
        end
            
                
        %% Solving
        function d = solveLinearElasticity(obj, f, t, sides, g, dofs)
            [K, ~] = obj.assembleMatrices;
            F = obj.assembleForce(f);
            for i=1:numel(sides)
                FN = obj.applyTraction(t{i}, sides(i));
                F = F+FN;
            end
            [d, ~, ~] = obj.dirichlet_linear_solve(K,F,g,dofs);
        end
        
        function d = solveDynamicProblem(obj, f, omega, dofs)
            [K, M] = obj.assembleMatrices;
            F = obj.assembleForce(f);
            C = obj.alpha*M +obj.beta*K;
            Kd = K +1i*omega*C -omega*omega*M;
            g = zeros(numel(dofs),1);
            [d, ~, ~] = obj.dirichlet_linear_solve(Kd, F, g, dofs);
        end
        
        function [d, omega, modal_solution] = naturalModes(obj, sides, num_modes)
            [K, M] = obj.assembleMatrices;
            b = obj.domain.extract_boundaries;
            id = obj.id_matrix;
            clamped_dofs = [];
            for i=1:numel(sides)
                side = sides(i);
                cpoints = b{side,4};
                dofs = id(cpoints,:);
                dofs = dofs(:);
                clamped_dofs = [clamped_dofs; dofs];
            end
            clamped_dofs = unique(clamped_dofs);
            unclamped_dofs = setdiff(1:length(K), clamped_dofs);
            clamped_dofs = sort(clamped_dofs(:), 'asc');
            K(clamped_dofs,:) = [];
            K(:, clamped_dofs) = [];
            M(clamped_dofs,:) = [];
            M(:, clamped_dofs) = [];
            [V, W] = eigs(K, M, num_modes,'sm');
            omega = sqrt(diag(W));
            
            [~, sz2] = size(V);
            d = zeros(numel(clamped_dofs) +numel(unclamped_dofs),sz2);
            d(unclamped_dofs,:) = V;
            modal_solution = cell(num_modes,1);
            for i=1:num_modes
                modal_solution{i} = Solution(obj, d(:,i));
            end
        end
        
        function FRFs = pointFRF(obj, f, omega_range, clamped_dofs, points)
            [K, M] = obj.assembleMatrices;
            F = obj.assembleForce(f);
            C = obj.alpha*M +obj.beta*K;
            FRFs = zeros(numel(points),numel(omega_range));
            for j=1:numel(omega_range)
                omega = omega_range(j);
                Kd = K +1i*omega*C -omega*omega*M;
                g = zeros(numel(clamped_dofs),1);
                [d, ~, ~] = obj.dirichlet_linear_solve(Kd, F, g, clamped_dofs);
                FRFs(:,j) = d(points);
            end
        end
        
        %% Utility Functions
    end
end