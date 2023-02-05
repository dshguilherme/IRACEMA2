classdef CahnHilliard < Assembler & handle
    
    properties
        
        lambda;
        theta;
        diffusivity;
        concentration;
        dt;
        sol;
    end
    
    methods
  
        function obj = CahnHilliard(lambda,theta,diffusivity,initial_c,domain)
            obj@Assembler("gauss", 1, domain);
            obj.lambda = lambda;
            obj.theta = theta;
            obj.diffusivity = diffusivity;
            obj.concentration = initial_c;
            obj.sol = Solution(obj,initial_c);
        end
        
        function c = elementConcentration(obj, elm, element)
            ind = elm(:,element);
            c = obj.concentration(ind);
        end
        
        function [mu, dmu, d2mu] = chemicalPotential(obj, c)
            mu = 0.5/obj.theta.*log(c./(1-c))+1-2.*c;
            mu = mu*3*obj.diffusivity;
            
            dmu = 0.5./obj.theta.*1./(c.*(1-c))-2;
            dmu = dmu*3*obj.diffusivity;
            
            d2mu=0.5./obj.theta.*(2.*c-1)/(c.*c.*(1-c).*(1-c));
            d2mu=d2mu*3*obj.diffusivity;
        end
        
        function [M, dM, d2M] = mobility(obj, c)
            M=c.*(1-c);
            dM=1-2.*c;
            d2M=repmat(-2,size(c));
        end
        
        function [K, F] = assembleSystem(obj)
            dt = obj.dt;
            d = obj.dimensions;
            [global_basis_index, element_local_mapping, element_ranges] = ...
              GetConnectivityArrays(obj.domain);
            [~, lm] = BuildGlobalLocalMatrices(element_local_mapping, d);

            [qp, qw] = obj.quad_rule;
            n_quad = length(qw);
            D = obj.diffusivity;
            lmbda = obj.lambda;
            ndof = max(max(element_local_mapping))*d;
            [nel_dof, nel] = size(element_local_mapping);
            K = zeros(ndof);
            F = zeros(ndof,1);
            for e=1:nel
                K_e = zeros(nel_dof);
                F_e = zeros(1,nel_dof);
                for n=1:n_quad
                    q = qp(n,:);
                    [R, dR, d2R, J, c, gradc, nablac] = CHShape(obj, ...
                        q, global_basis_index, element_local_mapping, element_ranges, e);
                    Jmod = abs(J*qw(n));
                    M_k = D*R.*(1-R);
                    M_f = D*c.*(1-c);
                    gradM_k = D*(1-2*R).*dR;
                    gradM_f = D*(1-2*c).*gradc;
                    mu_k = (0.5/obj.theta)*log(R./(1-R)) +1 -2*R;
                    mu_f = (0.5/obj.theta)*log(c./(1-c)) +1 -2*c;
                    gradmu_k = (0.5/obj.theta).*dR./(R -R.^2) -2*dR;
                    gradmu_f = (0.5/obj.theta).*gradc./(c -c.^2) -2*gradc;
                    
                    K_e = K_e +Jmod*(R*R') +0.5*Jmod*dt*(dR*((M_k.*gradmu_k)' +lmbda*(gradM_k.*d2R)') +lmbda*(d2R*d2R'));
                    F_e(:) = F_e(:) +Jmod*(R*c) +0.5*Jmod*dt*(dR*((M_f*gradmu_f)' +(gradM_f*lmbda*nablac)') +d2R*lmbda*nablac);
                end
                idx = lm(:,e)';
                K(idx,idx) = K(idx,idx) +K_e;
                F(idx) = F(idx) +F_e(:);
            end
            K = sparse(K);
            F = sparse(F);
        end
        
        function [ms, K_hat, F_hat] = imposeFlowBC(obj, K, F)
            bdries = obj.domain.extract_boundaries;
            Points = bdries(:,4);
            Points = cell2mat(Points);
            % Find corner points
            corners = [Points(1,1); Points(2,1); Points(1,end); Points(2,end)];
            jump_x = Points(3,2) - Points(3,1);
            jump_y = Points(1,2) - Points(1,1);
            x_up_points = setdiff(Points(1,:), corners)';
            x_down_points = setdiff(Points(2,:), corners)';
            y_up_points = setdiff(Points(3,:), corners)';
            y_down_points = setdiff(Points(4,:), corners)';
            
            msx = [x_up_points, x_up_points+jump_x; x_down_points, x_down_points-jump_x];
            msy = [y_up_points, y_up_points+jump_y; y_down_points, y_down_points-jump_y];
            
            ms = [msx; msy];
            corner_cases = corners +[jump_x, jump_y; -jump_x, jump_y; jump_x, -jump_y; -jump_x, -jump_y];
            
            idx = arrayfun(@(x) find(ms(:,1) == x), corner_cases(:));
            ms(idx,:) = [];
            msc = repmat(corners,[2 1]);
            msc = [msc, corner_cases(:)];
            ms = sort([msc; ms]);
            [global_basis_index, element_local_mapping, element_ranges] = ...
              GetConnectivityArrays(obj.domain);
            T = eye(size(K));
            for i=1:length(ms)
                T(ms(i,:), ms(i,:)) = [1 0; 1 0];
            end
            K_hat = T'*K*T;
            F_hat = T'*F;
            K_hat(ms(:,2),:) = [];
            K_hat(:,ms(:,2)) = [];
            F_hat(ms(:,2)) = [];
        end
        
        function [eB, eI] = computeEnergies(obj)
             d = obj.dimensions;
            [global_basis_index, element_local_mapping, element_ranges] = ...
              GetConnectivityArrays(obj.domain);
            [~, lm] = BuildGlobalLocalMatrices(element_local_mapping, d);

            [qp, qw] = obj.quad_rule;
            n_quad = length(qw);
            lmbda = obj.lambda;
            [~, nel] = size(element_local_mapping);
            eB = 0;
            eI = 0;
            for e=1:nel
                for n=1:n_quad
                        q = qp(n,:);
                        [~, ~, ~, J, c, gradc, ~] = CHShape(obj, ...
                            q, global_basis_index, element_local_mapping, element_ranges, e);
                        Jmod = abs(J*qw(n));
                        eB = eB+ Jmod*((0.5/obj.theta)*(c*log(c/(1-c)) +log(1-c)) -c^2 +c);
                        eI = eI+ Jmod*0.5*lmbda*norm(gradc)*norm(gradc);
                end
            end
           end
        
        function all_c = solveSystem(obj, initial_dt, max_steps)
            obj.dt = initial_dt;
            i = 0;
            bad_ending = 0;
            all_c = zeros(length(obj.concentration),max_steps);
            while i < max_steps
                [K,F] = obj.assembleSystem;
                [ms, K_hat, F_hat] = obj.imposeFlowBC(K, F);
%                 [L,U] = ilu(K_hat,struct('type','ilutp','droptol',1e-6));
                [c,flag,relres,iter] = gmres(K_hat,F_hat,[],1e-6,40); %,L,U);
                if flag
                    obj.dt = 0.5*obj.dt;
                    bad_ending = bad_ending +1;
                    if bad_ending > 10
                        break
                    end
                elseif iter <= 30
                    obj.dt = 1.1*obj.dt;
                    obj.concentration = c;
                    i = i+1;
                    fprintf('Step: %i | dt: %.2e | Relative residual: %.2e \n', i,obj.dt,relres);
                    bad_ending = 0;
                    all_c(:,i) = c;
                elseif iter > 30
                    obj.concentration = c;
                    i = i+1;
                    fprintf('Step: %i | dt: %.2e | Relative residual: %.2e \n', i,obj.dt,relres);
                    bad_ending = 0;
                    all_c(:,i) = c;
                end
            end
        end
        
%         function residual = assembleResidual(obj)
%             d = obj.dimensions;
%               [global_basis_index, element_local_mapping, element_ranges] = ...
%                 GetConnectivityArrays(obj.domain);
%             [~, lm] = BuildGlobalLocalMatrices(element_local_mapping, d);
%             
%             [qp, qw] = obj.quad_rule;
%             n_quad = length(qw);
%             
%             ndof = max(max(element_local_mapping))*d;
%             [nel_dof, nel] = size(element_local_mapping);
%             alpha = 1/(obj.lambda);
%             RR = zeros(ndof,1);
%             for e=1:nel
%                 R_e = zeros(1,nel_dof);
%                 for n=1:n_quad
%                     q = qp(n,:);
%                     [R, dR, d2R, J, c, gradc, nablac] = CHShape(obj, ...
%                         q, global_basis_index, element_local_mapping, element_ranges, e);
%                     Jmod = abs(J*qw(n));
%                     R_e(:) = R_e(:)+ R*c_dot +dR*(alpha/3)*(gradc -6*c*(1-c)*gradc);
%                     R_e(:) = R_e(:)+ dR*(1-2*c)*nablac*gradc +d2R*c*(1-c)*nablac;
%                     R_e(:) = R_e(:)*Jmod;
%                 end
%                 idx = lm(:,e)';
%                 RR(idx) = RR(idx)+R_e(:);
%             end
%             residual = sparse(RR);
%         end
        
%         function tangents = assembleTangentMatrix(obj)
%             d = obj.dimensions;
%               [global_basis_index, element_local_mapping, element_ranges] = ...
%                 GetConnectivityArrays(obj.domain);
%             [~, lm] = BuildGlobalLocalMatrices(element_local_mapping, d);
%             
%             [qp, qw] = obj.quad_rule;
%             n_quad = length(qw);
%             
%             ndof = max(max(element_local_mapping))*d;
%             [nel_dof, nel] = size(element_local_mapping);
%             alpha = 1/(obj.lambda);
%             M = zeros(ndof);
%             K = zeros(ndof);
%             for e=1:nel
%                 K_m = zeros(nel_dof);
%                 K_k = zeros(nel_dof);
%                 for n=1:n_quad
%                     q = qp(n,:);
%                     [R, dR, d2R, J, c, ~, gradc, nablac] = CHShape(obj, ...
%                         q, global_basis_index, element_local_mapping, element_ranges, e);
%                     Jmod = abs(J*qw(n));
%                     K_m = K_m +Jmod*(R*R');
%                     K_k = K_k +dR*(alpha/3)*dR' -2*alpha*dR*(c*dR' +gradc*R' -(c^2)*dR' -2*c*gradc*R');
%                     K_k = K_k +dR*(gradc*d2R' +nablac*dR' -2*(c*gradc*d2R' +c*nablac*dR' +gradc*nablac*R'));
%                     K_k = K_k +d2R*(c*d2R' +nablac*R' -(c^2)*d2R' -2*c*nablac*R');
%                     K_k = K_k +dR*(alpha/3)*dR' -2*alpha*dR*((c.*dR)' +(gradc.*dR)' ...
%                         -((c.^2).*dR)' -2*(c.*(gradc.*dR))');
%                     K_k = K_k +dR*((gradc.*d2R)' + (gradc.*dR)') ...
%                         -2*dR*((c.*(gradc.*d2R))' +(c.*nablac.*dR)' +(gradc.*(nablac.*R))');
%                     K_k = K_k +d2R*((c.*d2R)' +(nablac.*R)' -((c.^2).*d2R)' -(2*c.*nablac.*R)');
%                     K_k = Jmod*K_k;
%                 end
%                 idx = lm(:,e)';
%                 M(idx,idx) = M(idx,idx)+K_m;
%                 K(idx,idx) = K(idx,idx)+K_k;
%             end
%             M = sparse(M);
%             K = sparse(K);
%             tangents = {M, K};
%         end

%         function c = solve_system(obj,alpha_m,alpha_f, gamma,dt)
%             obj.delta_t = dt;
%             obj.alpha_m = alpha_m;
%             obj.alpha_f = alpha_f;
%             obj.gamma = gamma;
% 
% %           predictor stage
%             Y_previous = obj.concentration;
%             Y_next = Y_previous;
%             Ydot_previous = obj.cdot;
%             Ydot_next = ((gamma-1)/gamma)*Ydot_previous;
%             
%             cond = 1;
%             i = 0;
%             initial_residual = obj.assembleResidual;
%             while (cond > 0)
%                 i = i+1;
%                 Y_between = Y_previous +alpha_f*(Y_next -Y_previous);
%                 Ydot_between = Ydot_previous +alpha_m*(Ydot_next -Ydot_previous);
%                 obj.concentration = Y_between;
%                 obj.cdot = Ydot_between;
%                 residual = obj.assembleResidual;
%                 tangents = obj.assembleTangentMatrix;
%                 M = tangents{1};
%                 K = tangents{2};
%                 K1 = alpha_m*M +alpha_f*gamma*obj.delta_t*K;
%                 deltaYdot = -residual\K1;
%                 Y_previous = Y_next;
%                 Ydot_previous = Ydot_next;
%                 Y_next = Y_next +gamma*obj.delta_t*(deltaYdot');
%                 Ydot_next = Ydot_next +deltaYdot';
%                 norm(residual)
%                 if norm(residual) <= (1e-4)*norm(initial_residual)
%                     cond = -1;
%                 end
%             end
%         end
        
    end
end