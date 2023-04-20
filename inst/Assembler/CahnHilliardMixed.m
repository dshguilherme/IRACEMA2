classdef CahnHilliardMixed < Assembler & handle
    properties
        lambda; % Interface thickness parameter
        mobility; % Diffusibility parameter
        theta; % Generalized theta-method integration 1 = BackwardsEuler, 0.5 = Crank-Nicolson
        c0; % Previous step concentration
        c1; % Next step concentration
        mu0; % Previous step chemical potential
        mu1; % Next step chemical potential
        dt; % Time discretization
        fc; % Free Energy
        dfdc; % Derivative of the Free Energy
        d2fdc; % Second derivative of the Free Energy
        timetable; % Data array for time integration
        solution_array; % Solution storage
        max_timesteps; % Maximum number of timesteps
        t_max; % Maximum time
        res_tol; % Acceptable tolerance for non-linear iteration
        newton_max_steps; % Maximum steps on newton iterations
        mode; % Native division or  preconditioned GMRES
        gmres_tol = 1e-5;
        gmres_maxit = 30;
    end
    
    methods
             
        %% Constructor for the class
        function obj = CahnHilliardMixed(lambda, mobility, domain, theta, dt, t_max, max_timesteps)
            obj@Assembler("gauss",1,domain);
            obj.lambda = lambda;
            obj.mobility = mobility;
            obj.theta = theta;
            obj.c0 = obj.initial_concentration;
            obj.solution_array = zeros(length(obj.c0), max_timesteps);
            obj.solution_array(:,1) = obj.c0;
            obj.dt = dt;
            obj.timetable = zeros(max_timesteps,6); %time, dt, eTotal, eBulk, eInterface, residual norm
            obj.fc = @(c) (c^2)*((1-c)^2);
            obj.dfdc = @(c) (2*c)*(2*(c^2) -3*c +1);
            obj.d2fdc = @(c) 2*(6*c^2 -6*c -1);
            obj.mu0 = zeros(size(obj.c0));
            obj.max_timesteps = max_timesteps;
            obj.t_max = t_max;
            [eT, eB, eI] = obj.computeEnergies;
            obj.timetable(1,:) = [0, obj.dt, eT, eB, eI, 1];
            obj.res_tol = 1e-4; % default residual tolerance
            obj.newton_max_steps = 30; % default maximum number of newton steps
            obj.mode = "Newton";
        end
        
        % Stochastic Initial conditions. Might want to change this for
        % non-academic applications
        function c0 = initial_concentration(obj)
            c0 = 0.63 +0.05*randn(numel(obj.domain.points),1);
        end
        
        % Projection of the initial conditions to mu
%         Maybe not necessary.
%         function mu0 = initial_potential(obj)
%             [gbi, elm, er, lm, qp, qw, n_quad, ndof, nel_dof, nel] = obj.preLoopParser;
%             K = zeros(ndof);
%             F = zeros(ndof,1);
%             for e=1:nel
%                 k = zeros(nel_dof);
%                 f = zeros(nel_dof,1);
%                 for n=1:n_quad
%                     q = qp(n,:);
%                     [R, dR, J] = FastShape(obj.domain, q, gbi, elm, er, e);
%                     Jmod = abs(J*qw(n));
%                     
%                     [c, gradc, ~, ~] = obj.localConcentrationInfo(R, dR, elm, e, obj.c0, obj.c0);
%                     k = k + Jmod*(R*R');
%                     f = f + Jmod*(obj.dfdc(c) +obj.lambda*(dR*gradc'));
%                 end
%                 idx = lm(:,e)';
%                 K(idx, idx) = K(idx, idx) + k;
%                 F(idx) = F(idx) +f;
%             end
%             mu0 = sparse(K)\sparse(F);
%         end
        
        %% Residual and tangent assembly
        function residual = assembleResidual(obj)
            [gbi, elm, er, lm, qp, qw, n_quad, ndof, nel_dof, nel] = obj.preLoopParser;
            residual = zeros(ndof*2,1);
            for e=1:nel
                r_ec = zeros(nel_dof,1);
                r_em = r_ec;
                for n=1:n_quad
                    q = qp(n,:);
                    [R, dR, J] = FastShape(obj.domain, q, gbi, elm, er, e);
                    Jmod = abs(J*qw(n));
                    
                    [c_0, ~, ~, gradmu0] = obj.localConcentrationInfo(R, dR, elm, e, obj.c0, obj.mu0);
                    [c_1, gradc1, mu_1, gradmu1] = obj.localConcentrationInfo(R, dR, elm, e, obj.c1, obj.mu1);
                    
                    grad_mumid = (1-obj.theta)*gradmu0 +obj.theta*gradmu1;
                    df = obj.dfdc(c_1);
                    M = obj.mobility;
                    
                    r_ec = r_ec + Jmod*(R*(c_1 - c_0) +obj.dt*M*(dR*grad_mumid'));
                    r_em = r_em +Jmod*(R*(mu_1 - df) -obj.lambda*(dR*gradc1'));
                end
                idx = lm(:,e)';
                residual(idx) = residual(idx) +r_ec;
                residual(ndof+idx) = residual(ndof+idx) +r_em;
            end
            residual = sparse(residual);
        end
        
        function tangent = assembleTangent(obj)
            [gbi, elm, er, lm, qp, qw, n_quad, ndof, nel_dof, nel] = obj.preLoopParser;
            Kc = zeros(ndof);
            Km = Kc;
            Kcm = Kc;
            Kmc = Kc;
            for e=1:nel
                k_cc = zeros(nel_dof);
                k_cm = k_cc;
                k_mc = k_cc;
                k_mm = k_cc;
                for n=1:n_quad
                    q = qp(n,:);
                    [R, dR, J] = FastShape(obj.domain, q, gbi, elm, er, e);
                    Jmod = abs(J*qw(n));
                    M = obj.mobility;
                    
                    [c_1, ~, ~, ~] = obj.localConcentrationInfo(R, dR, elm, e, obj.c1, obj.mu1);
%                     df2 = obj.d2fdc(c_1);
                    k_cc = k_cc + Jmod*(R*R');
                    k_cm = k_cm +Jmod*(M*(dR*dR'))*obj.dt;
%                     k_mc = k_mc +Jmod*(-df2*(R*R') -obj.lambda*(dR*dR'));
                    k_mc = k_mc +Jmod*(-obj.lambda*(dR*dR'));
                    k_mm = k_mm +Jmod*(R*R');
                end
                idx = lm(:,e)';
                Kc(idx,idx) = Kc(idx,idx) + k_cc;
                Kcm(idx,idx) = Kcm(idx,idx) +k_cm;
                Kmc(idx,idx) = Kmc(idx,idx) +k_mc;
                Km(idx,idx) = Km(idx,idx) + k_mm;
            end
            tangent = sparse([Kc, Kcm; Kmc, Km]);
        end
        
        %% Time Integration
        function [converged, n_steps, res_norm] = timeStep(obj, max_corrections)
            ndof = length(obj.c0);
            converged = 0;
            % Predictor Step
            obj.c1 = obj.c0;
            obj.mu1 = obj.mu0;
            % Multicorrector Step
            for i=1:max_corrections+1
                residual = obj.assembleResidual;
                res_norm = norm(residual);
                if res_norm < obj.res_tol
                    converged = 1;
                    break
                end
                tangent = obj.assembleTangent;
                switch obj.mode
                    case "Newton"
                        u = tangent\(-residual);
                    case "GMRES"
                        [L, U] = ilu(tangent,struct('type','ilutp','droptol',1e-6));
                        [u, flag, ~, ~, ~] = gmres(tangent,-residual,[],obj.gmres_tol,obj.gmres_maxit);
                        if flag
                            converged = 0;
                            break
                        end
                end
                obj.c1 = obj.c1 +u(1:ndof);
                obj.mu1 = obj.mu1 +u(ndof+1:end);        
            end
            n_steps = i;
        end
        
        function timeLoop(obj)
            t = 0;
            steps = 1;
            while (t < obj.t_max) && (steps < obj.max_timesteps)
                t = t+obj.dt;
                [converged, n_steps, res_norm] = obj.timeStep(obj.newton_max_steps);
                if converged
                    obj.c0 = obj.c1;
                    [eT, eB, eI] = obj.computeEnergies;
                    obj.timetable(steps,:) = [t, obj.dt, eT, eB, eI, res_norm];
                    obj.solution_array(:,steps) = obj.c1;
                    fprintf("Step: %i | Time: %.2e | dt: %.2e | Residual Norm: %.2e \n", steps, t, obj.dt , res_norm)
                    if n_steps < 7
                        obj.dt = 1.1*obj.dt;
                    end
                    steps = steps+1;
                else
                    obj.c1 = obj.c0;
                    t = t-obj.dt;
                    obj.dt = 0.5*obj.dt;
                end
            end
            obj.timetable = obj.timetable(1:steps-1,:);
            obj.solution_array = obj.solution_array(:,1:steps-1);
        end
           %% Aux Functions
           function plotEvolution(obj)
                t_array = obj.timetable(:,1);
                dt_array = obj.timetable(:,2);
                eT = obj.timetable(:,3);
                eB = obj.timetable(:,4);
                eI = obj.timetable(:,5);
                figure(3)
                loglog(t_array,dt_array)
                set(gca, 'FontSize',18)
                ylabel('dt [s]', 'FontWeight', 'bold', 'FontSize', 23)
                xlabel('time [s]', 'FontWeight', 'bold', 'FontSize', 23)
                title('Time Adaptativity', 'FontWeight', 'bold', 'FontSize', 23)
                grid on
                
                figure(4)
                hold on
                plot(t_array,eT, 'LineWidth', 3)
                plot(t_array,eB, 'LineWidth', 3)
                plot(t_array,eI, 'LineWidth', 3)
                legend('Total Energy', 'Bulk Energy', 'Interface Energy')
                set(gca, 'FontSize', 18)
                ylabel('Energy [J]', 'FontWeight', 'bold', 'FontSize', 23)
                xlabel('time [s]', 'FontWeight', 'bold', 'FontSize', 23)
                title('Free Energy Functional', 'FontWeight', 'bold', 'FontSize', 23)
                grid on
           end
           
            function [c, gradc, mu, grad_mu] = localConcentrationInfo(obj, R, dR, elm, e, c_vec, mu_vec)
                ind = elm(:,e);
                c = dot(R, c_vec(ind));
                mu = dot(R, mu_vec(ind));

                c_x = dot(c_vec(ind), dR(:,1));
                c_y = dot(c_vec(ind), dR(:,2));
                c_z = dot(c_vec(ind), dR(:,3));
                gradc = [c_x, c_y, c_z];

                mu_x = dot(mu_vec(ind), dR(:,1));
                mu_y = dot(mu_vec(ind), dR(:,2));
                mu_z = dot(mu_vec(ind), dR(:,3));
                grad_mu = [mu_x, mu_y, mu_z];
            end
    
        function [eT, eB, eI] = computeEnergies(obj)
            [gbi, elm, er, ~, qp, qw, n_quad, ~, ~, nel] = obj.preLoopParser;
            
            eB = 0;
            eI = 0;
            for e=1:nel
                for n=1:n_quad
                    q = qp(n,:);
                    [R, dR, J] = FastShape(obj.domain, q, gbi, elm, er, e);
                    Jmod = abs(J*qw(n));
                    
                    [c, gradc, ~, ~] = obj.localConcentrationInfo(R, dR, elm, e, obj.c0, obj.mu0);
                    
                    eB = eB+ Jmod*obj.fc(c);
                    eI = eI+ Jmod*obj.lambda*dot(gradc, gradc);
                end
            end
            eT = eB +eI;
        end
        
    end
            
end