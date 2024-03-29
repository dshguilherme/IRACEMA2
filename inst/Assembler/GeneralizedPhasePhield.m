classdef GeneralizedPhasePhield < CahnHilliardMixed & handle
   
    properties
        eta; % Field coupling constant
        cf_asb; % Assembler of the coupled field
        d0; % Previous solution of the coupled field
        d1; % Next solution of the coupled field
        trac; % Traction vector on the coupled field
        sides; % Sides of the traction
        clamped_dofs; % dofs that are clamped
        force_function; % Function of (x) of the body force
        frequency; % Frequency of the problem
    end
    
    methods
        
        function obj = GeneralizedPhasePhield(coupled_field_assembler, trac, ...
                sides, clamped_dofs, force_function, freq, eta, lambda, mobility, domain, theta, dt, t_max, max_timesteps)
            obj@CahnHilliardMixed(lambda, mobility, domain, theta, dt, t_max, max_timesteps);
            obj.cf_asb = coupled_field_assembler;
            obj.eta = eta;
            obj.trac = trac;
            obj.sides = sides;
            obj.clamped_dofs = clamped_dofs;
            obj.force_function = force_function;
            obj.frequency = freq;
            tmp = obj.timetable(1,:);
            obj.timetable = zeros(max_timesteps,8); %time, dt, eT, eB, eI, eP, eK, res_norm
            obj.timetable(1,1:6) = tmp;
        end
        
        %% Assembly
        
        function residual = assembleStaggeredResidual(obj, option)
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
                    [u, grad_u] = obj.cf_asb.localDisplacementInfo(R, dR, elm, e, obj.d0);
                    D = obj.cf_asb.D0;
                    switch obj.cf_asb.dimensions
                        case 2
                            epsilon = [diag(grad_u); grad_u(1,2) + grad_u(2,1)];
                        case 3
                            epsilon = [diag(grad_u); grad_u(2,3) + grad_u(3,2); ...
                                grad_u(3,1) + grad_u(1,3); grad_u(2,1) +grad_u(1,2)];
                    end
                    ds = real(dot(3*c_1*c_1*D*epsilon, epsilon));
                    dk = 3*c_1*c_1*obj.cf_asb.rho*dot(u,u);
                    
                    r_ec = r_ec + Jmod*(R*(c_1 - c_0) +obj.dt*M*(dR*grad_mumid'));
                    switch option
                        case "elastic"
                            r_em = r_em +Jmod*(R*(mu_1 - df +ds) -obj.lambda*(dR*gradc1'));
                        case "dynamic"
                            r_em = r_em +Jmod*(R*(mu_1 - df +obj.eta*(ds -dk)) -obj.lambda*(dR*gradc1'));
                    end
                end
                idx = lm(:,e)';
                residual(idx) = residual(idx) +r_ec;
                residual(ndof+idx) = residual(ndof+idx) +r_em;
            end
            residual = sparse(residual);            
        end
        
        function tangent = assembleStaggeredTangent(obj, option)
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
                    df2 = obj.d2fdc(c_1);
                    [u, grad_u] = obj.cf_asb.localDisplacementInfo(R, dR, elm, e, obj.d0);
                    D = obj.cf_asb.D0;
                    switch obj.cf_asb.dimensions
                        case 2
                            epsilon = [diag(grad_u); grad_u(1,2) + grad_u(2,1)];
                        case 3
                            epsilon = [diag(grad_u); grad_u(2,3) + grad_u(3,2); ...
                                grad_u(3,1) + grad_u(1,3); grad_u(2,1) +grad_u(1,2)];
                    end
                    ds = real(dot(6*c_1*D*epsilon, epsilon));
                    dk = 6*c_1*obj.cf_asb.rho*dot(u,u);
                    
                    k_cc = k_cc + Jmod*(R*R');
                    k_cm = k_cm +Jmod*(M*(dR*dR'))*obj.dt;
                    switch option
                        case "elastic"
                            k_mc = k_mc +Jmod*((obj.eta*ds -df2)*(R*R') -obj.lambda*(dR*dR'));
                        case "dynamic"
                            k_mc = k_mc +Jmod*((obj.eta*(ds-dk) -df2)*(R*R') -obj.lambda*(dR*dR'));
                    end
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
        
        function d = solveCoupledField(obj, option)
            g = zeros(size(obj.clamped_dofs));
            switch option
                case "elastic"
                    d = obj.cf_asb.solveConcentratedElasticity(obj, obj.force_function, ...
                        obj.trac, obj.sides, g, obj.clamped_dofs);
                case "dynamic"
                    d = obj.cf_asb.solveConcentratedDynamic(obj, obj.force_function, ...
                        obj.frequency, obj.clamped_dofs);
            end
        end
        
        function [residual, tangent] = assembleMonolythicSystem(obj)
        end
        
        %% Time Integration
        function staggeredTimeLoop(obj, option)
            obj.d0 = obj.solveCoupledField(option);
            t = 0;
            steps = 1;
            while (t < obj.t_max) && (steps < obj.max_timesteps)
                t = t+obj.dt;
                [converged, n_steps, res_norm] = obj.staggeredTimeStep(obj.newton_max_steps, option);
                if converged
                    obj.c0 = obj.c1;
                    obj.solution_array(:,steps) = obj.c1;
                    [~, eB, eI] = obj.computeEnergies;
                    obj.d0 = obj.solveCoupledField(option);
                    [eP, eK] = obj.cf_asb.computeEnergies(obj.frequency, obj.d0);
                    switch option
                        case "elastic"
                            eT = eB+eI+eP;
                            obj.timetable(steps,:) = [t, obj.dt, eT, eB, eI, eP, 0, res_norm];
                        case "dynamic"
                            eT = eB+eI+eP+eK;
                            obj.timetable(steps,:) = [t, obj.dt, eT, eB, eI, eP, eK, res_norm];
                    end
                    fprintf("Step: %i | Time: %.2e | dt: %.2e | Residual Norm: %.2e \n", steps, t, obj.dt , res_norm)
                    if n_steps < 5
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
        
        function [converged, n_steps, res_norm] = staggeredTimeStep(obj, max_corrections, option)
            ndof = length(obj.c0);
            converged = 0;
            % Predictor Step
            obj.c1 = obj.c0;
            obj.mu1 = obj.mu0;
            % Multicorrector Step
            for i=1:max_corrections+1
                residual = obj.assembleStaggeredResidual(option);
                res_norm = norm(residual);
                if res_norm < obj.res_tol
                    converged = 1;
                    break
                end
                tangent = obj.assembleStaggeredTangent(option);
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
        
        function monolythicTimeLoop(obj)
        end
        
        %% Auxiliary Functions
        function plotEvolution(obj, option)
                t_array = obj.timetable(:,1);
                dt_array = obj.timetable(:,2);
                eT = obj.timetable(:,3);
                eB = obj.timetable(:,4);
                eI = obj.timetable(:,5);
                eP = obj.timetable(:,6);
                if option == "dynamic"
                    eK = obj.timetable(:,7);
                end
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
                plot(t_array,eP, 'LineWidth', 3)
                if option == "dynamic"
                    plot(t_array,eK, 'LineWidth', 3)
                    legend('Total Energy', 'Bulk Energy', 'Interface Energy', 'Elastic Energy', 'Kinetic Energy')
                else
                    legend('Total Energy', 'Bulk Energy', 'Interface Energy', 'Elastic Energy')
                end
                set(gca, 'FontSize', 18)
                ylabel('Energy [J]', 'FontWeight', 'bold', 'FontSize', 23)
                xlabel('time [s]', 'FontWeight', 'bold', 'FontSize', 23)
                title('Free Energy Functional', 'FontWeight', 'bold', 'FontSize', 23)
                grid on
           end
        
        
        
    end

end