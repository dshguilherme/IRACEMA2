classdef CahnHilliard < Assembler & handle
    
    properties
        lambda;
        theta;
        diffusivity;
        concentration;
        cdot;
        dt;
        mixed;
        mu;
        timetable;
        penalty;
    end
    
    methods
  
        function obj = CahnHilliard(lambda,theta,diffusivity,initial_c,domain,mixed, max_timesteps)
            obj@Assembler("gauss", 1+mixed, domain);
            obj.lambda = lambda;
            obj.theta = theta;
            obj.diffusivity = diffusivity;
            obj.cdot = zeros(size(initial_c));
            obj.concentration = zeros(length(initial_c), max_timesteps);
            obj.concentration(:,1) = initial_c;
            obj.cdot = zeros(size(obj.concentration));
            obj.mixed = mixed;
            
            obj.dt = 5e-8;
            obj.timetable = zeros(max_timesteps,5); % time, dt, eTotal, eBulk, eInterface
            [eT, eB, eI] = obj.computeEnergies(1);
            obj.timetable(1,:) = [0, obj.dt, eT, eB, eI];
        end
        
        function [cn1, cdotn1] = alphaStep(obj, cn, cdotn, alpha_m, alpha_f, gamma)
%             Predictor Stage
            cn1 = cn;
            cdotn1 = ((gamma-1)/gamma)*cdotn;
%             Multicorrector Stage
            for i=1:100
                cdot_alpha = cdotn1 +alpha_m*(cdotn1 - cdotn);
                c_alpha = cn1 +alpha_f*(cn1-cn);
                
                [tangent, residual] = obj.assemblePrimalSystem(c_alpha,cdot_alpha,alpha_f, alpha_m, gamma);
                if norm(residual) < 1e-5
                    break
                end
                
                [L,U] = ilu(tangent,struct('type','ilutp','droptol',1e-6));
                [deltacdot, ~, ~, ~, ~] = gmres(tangent,-residual,[],1e-4,100, L, U);
                
                cdotn1 = cdotn1 +deltacdot;
                cn1 = cn1 +gamma*obj.dt*deltacdot;
                if i == 100
                    warning('Maximum iterations reached in alphaStep without reaching the desired tolerance.');
                end
            end
        end
        
        
        function c = backwardsEuler(obj, c, cdot_initial)
            if nargin < 3
                cdot_initial = zeros(size(c));
            end
            [c, ~] = obj.alphaStep(c, cdot_initial, 1, 1, 1); 
        end

        
        function c = elementConcentration(obj, elm, element)
            ind = elm(:,element);
            c = obj.concentration(ind);
        end
        
        function [mu, dmu, d2mu] = chemicalPotential(obj, c)
            mu = 0.5/obj.theta.*log(c./(1-c))+1-2.*c;
            mu = mu/obj.lambda;
            
            dmu = 0.5./obj.theta.*1./(c.*(1-c))-2;
            dmu = dmu/obj.lambda;
            
            d2mu=0.5./obj.theta.*(2.*c-1)/(c.*c.*(1-c).*(1-c));
            d2mu=d2mu/obj.lambda;
        end
        
        function [M, dM, d2M] = mobility(obj, c)
            M=c.*(1-c);
            dM=1-2.*c;
            d2M=repmat(-2,size(c));
        end
        
%         function [tangent, residual] = assembleMixedSystem(obj, cmid, mumid)
%            d = obj.dimensions;
%            [global_basis_index, element_local_mapping, element_ranges] = ...
%                GetConnectivityArrays(obj.domain);
%            [~, lm] = BuildGlobalLocalMatrices(element_local_mapping, d);
%            
%            [qp, qw] = obj.quad_rule;
%            n_quad = length(qw);
%            
%            ndof = max(max(element_local_mapping))*d;
%            [nel_dof, nel] = size(element_local_mapping);
%            residual = zeros(ndof,1);
%            tangent = zeros(ndof);
%            for e=1:nel
%                r_e = zeros(nel_dof,d);
%                K_c = zeros(nel_dof);
%                K_mu = zeros(nel_dof);
%                K_cmu = K_c;
%                K_muc = K_c;
%                for n=1:n_quad
%                    q = qp(n,:);
%                    [R, dR, Jmod, c, gradc, mu, gradmu] = CHMixedShape(obj, cmid, mumid, ...
%                         q, global_basis_index, element_local_mapping, element_ranges, e);
%                    [M, dM, ~] = obj.mobility(c);
%                    [dfdc, d2fdc2, ~] = obj.chemicalPotential(c);
%                    
%                    r_e(:,1) = Jmod*qw(n)*(R*c +obj.dt*dR*(M*gradmu'));
%                    r_e(:,2) = Jmod*qw(n)*(R*(mu -dfdc) -dR*gradc');
%                    
%                    K_c = Jmod*qw(n)*(R*R' +dR*obj.dt*(dM*gradmu')*R');
%                    K_cmu = Jmod*qw(n)*(0.5*M*dR*dR');
%                    K_muc = -Jmod*qw(n)*0.5*(R*d2fdc2*R' -dR*dR');
%                    K_mu = Jmod*qw(n)*0.5*(R*R');
%                end
%                idx = lm(:,e)';
%                idx_c = idx(1:length(idx)/2);
%                idx_mu = setdiff(idx, idx_c);
%                residual(idx_c) = residual(idx_c) +r_e(:,1);
%                residual(idx_mu) = residual(idx_mu) +r_e(:,2);
%                
%                tangent(idx_c,idx_c) = tangent(idx_c,idx_c) +K_c;
%                tangent(idx_mu,idx_mu) = tangent(idx_mu,idx_mu) +K_mu;
%                tangent(idx_c,idx_mu) = tangent(idx_c,idx_mu) +K_cmu;
%                tangent(idx_mu,idx_c) = tangent(idx_mu,idx_c) +K_muc;
%            end
%            residual = sparse(residual);
%            tangent = sparse(tangent);
%         end
    function [tangent, residual] = assembleMixedSystem(obj, c, c0)
           d = obj.dimensions;
           [global_basis_index, element_local_mapping, element_ranges] = ...
               GetConnectivityArrays(obj.domain);
           [~, lm] = BuildGlobalLocalMatrices(element_local_mapping, d);
           
           [qp, qw] = obj.quad_rule;
           n_quad = length(qw);
           
           ndof = max(max(element_local_mapping))*d;
           [nel_dof, nel] = size(element_local_mapping);
           residual = zeros(ndof,1);
           tangent = zeros(ndof);
           for e=1:nel
               r_ec = zeros(nel_dof, 1);
               r_em = r_ec;
               K_cc = zeros(nel_dof);
               K_cm = K_cc;
               K_mm = K_cc;

               for n=1:n_quad
                   q = qp(n,:);
                   [R, dR, Jmod, c, gradc, mu, gradmu] = CHMixedShape(obj, cmid, mumid, ...
                        q, global_basis_index, element_local_mapping, element_ranges, e);
                   M = obj.diffusivity;
                   [dfdc, ~, ~] = obj.chemicalPotential(c);
                   
                   r_ec = R*c/obj.dt + M*dR*gradmu';
                   r_ec = Jmod*qw(n)*r_ec;
                   
                   r_em = R*(mu-dfdc) -obj.lambda*dR*gradc';
                   r_em = Jmod*qw(n)*r_em;
                   
                   K_cc = Jmod*qw(n)*(R*R');
                   K_cm = Jmod*qw(n)*(dR*dR'); % M*obj.dt*
                   K_mc = K_cm;
                   K_mm = Jmod*qw(n)*(R*R'); % -lambda*
               end
               idx = lm(:,e)';
               idx_c = idx(1:length(idx)/2);
               idx_mu = setdiff(idx, idx_c);
               residual(idx_c) = residual(idx_c) +r_ec;
               residual(idx_mu) = residual(idx_mu) +r_em;
               
               tangent(idx_c,idx_c) = tangent(idx_c,idx_c) +K_cc;
               tangent(idx_mu,idx_mu) = tangent(idx_mu,idx_mu) +K_mm;
               tangent(idx_c,idx_mu) = tangent(idx_c,idx_mu) +K_cm;
               tangent(idx_mu,idx_c) = tangent(idx_mu,idx_c) +K_mc;
           end
           residual = sparse(residual);
           tangent = sparse(tangent);
        end        
        function [tangent, residual] = assemblePrimalSystem(obj,cf,cdotm,alpha_f, alpha_m, gamma)
            dt = obj.dt;
            d = obj.dimensions;
           [global_basis_index, element_local_mapping, element_ranges] = ...
               GetConnectivityArrays(obj.domain);
           [~, lm] = BuildGlobalLocalMatrices(element_local_mapping, d);
           
           [qp, qw] = obj.quad_rule;
           n_quad = length(qw);
           
           ndof = max(max(element_local_mapping))*d;
           [nel_dof, nel] = size(element_local_mapping);
           residual = zeros(ndof,1);
           tangent = zeros(ndof);
           
           for e=1:nel
               r_e = zeros(nel_dof,d);
               K_e = zeros(nel_dof,d);
               for n=1:n_quad
                   q = qp(n,:);
                   [R, dR, d2R, Jmod, c, cdot, gradc, lapc] = CHShape(obj, cf, cdotm, ...
                       q, global_basis_index, element_local_mapping, element_ranges, e);
                   [~, dmu, d2mu] = obj.chemicalPotential(c);
                   [M, dM, d2M] = obj.mobility(c);
                   Jmod = Jmod*qw(n);
                   % Residual
                   r_e = r_e + Jmod*R*cdot;
                   r_e = r_e +Jmod*sum(dR.*(M*dmu*gradc +dM*gradc*lapc),2);
                   r_e = r_e +Jmod*d2R*(M*lapc);

                   % Tangent !
                   zeroth = R*R';
                   first = (dR*dR')*(M*dmu +lapc*dM);
                   second = (dR*gradc')*(M*d2mu + dmu*dM +d2M*lapc)*R';
                   third = (dR*gradc')*(dM)*d2R';
                   fourth = d2R*(lapc*dM)*R';
                   fifth = d2R*M*d2R';
                   K_e = K_e +Jmod*alpha_m*zeroth;
                   K_e = K_e +Jmod*alpha_f*gamma*obj.dt*(first+second+third+fourth+fifth);
               end
               idx = lm(:,e)';
               residual(idx) = residual(idx) +r_e(:);
               tangent(idx,idx) = tangent(idx,idx) +K_e;
           end
%            bdries = obj.domain.extract_boundaries;
%            normal = [-1 0 0; 1 0 0; 0 -1 0; 0 1 0;];
%            for i=1:length(bdries)
%                direction = setdiff([1 2], round(i/2));
%                elements = bdries{i,2};
%                [qp, qw] = gaussian_quadrature(bdries{i,1}.p);
%                n_quad = length(qw);
%                n_vec = normal(i,:);
%                for j=1:length(elements)
%                    e = elements(j);
%                    r_e = zeros(nel_dof, d);
%                    K_e = zeros(nel_dof);
%                    for n=1:n_quad
%                        q = zeros(1,2);
%                        q(round(i/2)) = (-1)^(mod(i,2));
%                        q(direction) = qp(n);
%                    [R, dR, d2R, Jmod, ~, ~, gradc, lapc] = CHShape(obj, cf, cdotm, ...
%                        q, global_basis_index, element_local_mapping, element_ranges, e);
%                    
%                    r_e = r_e +Jmod*qw(n)*(dR*n_vec')*(obj.penalty*dot(n_vec,gradc) -lapc);
%                    r_e = r_e -Jmod*qw(n)*d2R*(dot(n_vec,gradc));
%                    
%                    K_e = K_e +Jmod*qw(n)*((dR*n_vec')*(obj.penalty*n_vec*dR' -d2R') -d2R*(n_vec*dR'));
%                    end
%                    idx = lm(:,e)';
%                    residual(idx) = residual(idx) +r_e;
%                    tangent(idx,idx) = tangent(idx,idx) +K_e;
%                end
%            end
           residual = sparse(residual);
           tangent = sparse(tangent);
        end            
%         
%         function tangent = assemblePrimalTangent(obj, cf, cdotm)
%         end
%         
%         function stepTime(obj, i)
%             assert(i > 1, "You're trying to step the initial conditions. Please start loop at i=2");
%             if obj.mixed
%                 [c, mu] = obj.crankNicolsonStep;
%                 obj.concentration(:,i) = c(:);
%                 obj.cn = c(:);
%                 obj.mun = mu;
%  
%             else
%                 cond = 1;
%                 c = obj.cn;
%                 cdot = obj.cdotn;
%                 [eOld, ~, ~] = obj.computeEnergies;
%                 while cond > 0
%                     obj.cn = c;
%                     obj.cdotn = cdot;
%                     [cBE, ~] = obj.backwardsEuler;
%                     obj.cn = c;
%                     obj.cdotn = cdot;
%                     [cAlpha, cdotAlpha] = obj.generalizedAlphaStep;
%                     error = norm(cBE - cAlpha)/norm(cAlpha);
%                     tol = 1e-3;
%                     obj.cn = cAlpha;
%                     [eT, eB, eI] = obj.computeEnergies;
%                     if (error > tol)
%                         obj.dt = obj.dt*0.9*((tol/error)^0.5);
%                     else
%                         obj.dt = obj.dt*0.9*((tol/error)^0.5);
%                         obj.cn = cAlpha;
%                         obj.cdotn = cdotAlpha;
%                         cond = -1;
%                     end
%                 end
%                 obj.concentration(:,i) = cAlpha;
%             end
%             time = obj.timetable(i-1,1) +obj.dt;
%             obj.timetable(i,:) = [time, obj.dt, eT, eB, eI];
%         end
%         
%         function [c, cdot] = backwardsEuler(obj)
%             % predictor step
%             obj.cn1 = obj.cn;
%             obj.cdotn1 = zeros(size(obj.cn));
%             % multicorrector stage
%             for i=1:6
%                 [K, R] = obj.assemblePrimalSystem(obj.cn1, obj.cdotn1, 1, 1, 1);
% %                 [ms, K, R] = obj.imposeFlowBC(K,R);
%                 [L,U] = ilu(K,struct('type','ilutp','droptol',1e-6));
%                 [deltacdot,flag,relres,iter] = gmres(K,-R,[],1e-4,40,L,U);
% %                   deltacdot = (-R\K)';
% %                 deltacdot = ones(length(K)+length(ms),1);
% %                 deltacdot(ms(:,2)) = 0;
% %                 idx = find(abs(deltacdot -1) < eps);
% %                 deltacdot(idx) = dcdot;
% %                 deltacdot(ms(:,2)) = deltacdot(ms(:,1));
%                 obj.cdotn1 = obj.cdotn1 +deltacdot;
%                 obj.cn1 = obj.cn1 +obj.dt*deltacdot;
% %                 if relres <1e-4
% %                     break
% %                 end
%             end
%             c = obj.cn1;
%             cdot = obj.cdotn1;
%         end
%         
%         function [c, cdot] = generalizedAlphaStep(obj)
%             rho_inf = 0.5;
%             alpha_m = 0.5*((3-rho_inf)/(1+rho_inf));
%             alpha_f = 1/(1+rho_inf);
%             gamma = 0.5 +alpha_m -alpha_f;
%             % predictor step
%             obj.cn1 = obj.cn;
%             obj.cdotn1 = ((gamma-1)/gamma)*obj.cdotn;
%             
%             % multicorrector stage
%             for i=1:100
%                 cdotm = obj.cdotn +alpha_m*(obj.cdotn1 -obj.cdotn);
%                 cf = obj.cn +alpha_f*(obj.cn1 -obj.cn);
%                 [K,R] = obj.assemblePrimalSystem(cf, cdotm, alpha_f, alpha_m, gamma);
%                 if norm(R) < 1e-3
%                     break
%                 end
% %                 [ms, K, R] = obj.imposeFlowBC(K,R);
% %                 [L,U] = ilu(K,struct('type','ilutp','droptol',1e-6));
% %                 [deltacdot,flag,relres,iter] = gmres(K,-R,[],1e-4,40,L,U);
%                 deltacdot = K\(-R);
% %                   deltacdot = (-R\K)';
% %                 deltacdot = ones(length(K)+length(ms),1);
% %                 deltacdot(ms(:,2)) = 0;
% %                 idx = find(abs(deltacdot -1) < eps);
% %                 deltacdot(idx) = dcdot;
% %                 deltacdot(ms(:,2)) = deltacdot(ms(:,1));
%                 obj.cdotn1 = obj.cdotn1 +deltacdot;
%                 obj.cn1 = obj.cn1 +gamma*obj.dt*deltacdot;
% %                 if relres < 1e-4
% %                     break
% %                 end
%             end
%             c = obj.cn1;
%             cdot = obj.cdotn1;
%         end
%         
%         function [c, mu] = crankNicolsonStep(obj)
%             % Predictor stage
%             obj.cn1 = obj.cn;
%             obj.mun1 = obj.mun;
%             
%             % Multicorrector stage
%             for i=1:10
%                 mumid = 0.5*(obj.mun +obj.mun1);
%                 cmid = 0.5*(obj.cn +obj.cn1);
%                 [K, R] = obj.assembleMixedSystem(cmid, mumid);
% %                 [L,U] = ilu(K,struct('type','ilutp','droptol',1e-6));
% %                 [dc,flag,relres,iter] = gmres(K,-R,[],1e-4,40,L,U);
%                 dc = -R\K;
%                 dc = dc(:);
%                 deltac = dc(1:length(dc)/2);
%                 deltamu = dc(length(dc)/2+1:end);
%                 obj.mun1 = obj.mun1 +deltamu(:);
%                 obj.cn1 = obj.cn1 +deltac(:);
%             end
%             c = obj.cn1;
%             mu = obj.mun1;
%         end
%        
%         function [ms, K_hat, F_hat] = imposeFlowBC(obj, K, F)
%             bdries = obj.domain.extract_boundaries;
%             Points = bdries(:,4);
%             Points = cell2mat(Points);
%             % Find corner points
%             corners = [Points(1,1); Points(2,1); Points(1,end); Points(2,end)];
%             jump_x = Points(3,2) - Points(3,1);
%             jump_y = Points(1,2) - Points(1,1);
%             x_up_points = setdiff(Points(1,:), corners)';
%             x_down_points = setdiff(Points(2,:), corners)';
%             y_up_points = setdiff(Points(3,:), corners)';
%             y_down_points = setdiff(Points(4,:), corners)';
%             
%             msx = [x_up_points, x_up_points+jump_x; x_down_points, x_down_points-jump_x];
%             msy = [y_up_points, y_up_points+jump_y; y_down_points, y_down_points-jump_y];
%             
%             ms = [msx; msy];
%             corner_cases = corners +[jump_x, jump_y; -jump_x, jump_y; jump_x, -jump_y; -jump_x, -jump_y];
%             
%             idx = arrayfun(@(x) find(ms(:,1) == x), corner_cases(:));
%             ms(idx,:) = [];
%             msc = repmat(corners,[2 1]);
%             msc = [msc, corner_cases(:)];
%             ms = sort([msc; ms]);
% %             idx = find(ms(:,1) > ms(:,2));
% %             ms(idx,:) = fliplr(ms(idx,:));
%             [global_basis_index, element_local_mapping, element_ranges] = ...
%               GetConnectivityArrays(obj.domain);
%             T = eye(size(K));
%             for i=1:length(ms)
%                 T(ms(i,:), ms(i,:)) = [1 0; 1 0];
%             end
%             K_hat = T'*K*T;
%             F_hat = T'*F;
%             K_hat(ms(:,2),:) = [];
%             K_hat(:,ms(:,2)) = [];
%             F_hat(ms(:,2)) = [];
%             K_hat = sparse(K_hat);
%             F_hat = sparse(F_hat);
%         end
%         
        function [eT, eB, eI] = computeEnergies(obj, i)
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
            cf = obj.concentration(:,i);
            cdotm = obj.cdot(:,i);
            for e=1:nel
                for n=1:n_quad
                        q = qp(n,:);
                            [~, ~, ~, J, c, ~, gradc, ~] = CHShape(obj, cf, cdotm, ...
                                q, global_basis_index, element_local_mapping, element_ranges, e);
                            Jmod = abs(J*qw(n));
                            eB = eB+ Jmod*(c*log(c) +(1-c)*log(1-c) +2*obj.theta*c*(1-c));
                            eI = eI+ Jmod*0.5*obj.theta*lmbda*dot(gradc,gradc);
                            eT = eB+eI;
                end
            end
        end
%            
        function penalty = evalPenalty(obj)
             d = obj.dimensions;
            [global_basis_index, element_local_mapping, element_ranges] = ...
              GetConnectivityArrays(obj.domain);
            [~, lm] = BuildGlobalLocalMatrices(element_local_mapping, d);

            [qp, qw] = obj.quad_rule;
            n_quad = length(qw);
            ndof = max(max(element_local_mapping))*d;
           [nel_dof, nel] = size(element_local_mapping);
            A = zeros(ndof);
            B = zeros(ndof);
            for e=1:nel
                B_e = zeros(nel_dof);
                for n=1:n_quad
                        q = qp(n,:);
                        [~, ~, d2R, J, ~, ~, ~, ~] = CHShape(obj, obj.concentration(:,1), obj.cdot(:,1), ...
                            q, global_basis_index, element_local_mapping, element_ranges, e);
                        Jmod = abs(J*qw(n));
                        B_e = B_e + Jmod*d2R*d2R';
                end
                idx = lm(:,e)';
                B(idx,idx) = B(idx,idx) + B_e;
            end
           bdries = obj.domain.extract_boundaries;
           for i=1:length(bdries)
               direction = setdiff([1 2], round(i/2));
               elements = bdries{i,2};
               [qp, qw] = gaussian_quadrature(bdries{i,1}.p);
               n_quad = length(qw);
               for j=1:length(elements)
                   e = elements(j);
                   A_e = zeros(nel_dof);
                   for n=1:n_quad
                       q = zeros(1,2);
                       q(round(i/2)) = (-1)^(mod(i,2));
                       q(direction) = qp(n);
                   [~, ~, d2R, Jmod, ~, ~, ~, ~] = CHShape(obj, obj.concentration(:,1), obj.cdot(:,1), ...
                       q, global_basis_index, element_local_mapping, element_ranges, e);
                   A_e = A_e +Jmod*qw(n)*d2R*d2R';
                   end
                   idx = lm(:,e)';
                   A(idx,idx) = A(idx,idx) +A_e;
               end
           end
           d = eigs(A,B,1);
           penalty = 2.1*max(d);
        end
        
    
        
    end
end