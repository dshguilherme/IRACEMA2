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
        inter_geo; % Topological interpolation function domain
        D; % array of adimensionals for the problem
    end
    
    methods
        %% Initialization
        function obj = GeneralizedPhasePhield(coupled_field_assembler, trac, h0, ...
                sides, clamped_dofs, force_function, f0, freq, eta, lambda, ...
                mobility, domain, theta, dt, t_max, max_timesteps)
            obj@CahnHilliardMixed(lambda, mobility, domain, theta, dt, ...
                t_max, max_timesteps);
            obj.cf_asb = coupled_field_assembler;
            obj.eta = eta;
            obj.trac = trac;
            obj.sides = sides;
            obj.clamped_dofs = clamped_dofs;
            obj.force_function = force_function;
            obj.frequency = freq;
            tmp = obj.timetable(1,:);
            obj.timetable = zeros(max_timesteps,7); %time, dt, eT, eB, eI, eP, res_norm
            obj.timetable(1,1:6) = tmp;
            
            U = [0.055 0.055 0.055 0.055, 0.065, 0.075, .925, .935, .945, .945 , .945, .945];
            points = {[0 0 0 1]; [0 0 0 1]; [0 0 0 1]; [0 0 0 1]; [1 0 0 1]; [1 0 0 1]; [1 0 0 1]; [1 0 0 1]};
            p = 3;
            obj.inter_geo = Geometry(1, {U}, points, p);
            obj.c0 = obj.initial_concentration;
            obj.D = zeros(5,1);
            obj.D(4) = obj.cf_asb.young_modulus/h0;
            obj.D(5) = f0/h0;
            obj.D = obj.setAdimensionals(h0, f0);
        end
        
        function d0 = initial_displacement(obj)
            free_dofs = setdiff(1:length(obj.c0)*2,obj.clamped_dofs);
            [K, ~, F] = obj.assembleElasticitySystem;
            obj.d0 = zeros(length(obj.c0)*2,1);
            obj.d0(free_dofs) = K(free_dofs,free_dofs)\F(free_dofs);
            d0 = obj.d0;
        end
         
        function D = setAdimensionals(obj, h0, f0)
            [h, A] = obj.characteristic_length;
            L0 = A^(1/(obj.domain.rank));
            M0 = obj.mobility;
            lambda_ = obj.lambda*(h^2);
            T0 = (L0^4)/(lambda_*M0);
            E0 = obj.cf_asb.young_modulus;
            D2 = (L0^2)/(lambda_*h^2);
            D2_ = D2*lambda_;
            eta = obj.calculateAdimEta(D2_, E0);
            D3 = eta*obj.eta*D2;
            D4 = E0/h0;
            D5 = f0*L0/h0;
            D = [(L0^2); D2; D3; D4; D5];
        end
        
        function eta = calculateAdimEta(obj, D2, E0)
            d0 = obj.initial_displacement;
            [gbi, elm, er, lm, qp, qw, n_quad, ndof, nel_dof, nel] = obj.preLoopParser;
            for e=1:nel
                num = 0;
                den = 0;
                for n=1:n_quad
                    q = qp(n,:);
                    [R, dR, J] = FastShape(obj.domain, q, gbi, elm, er, e);
                    Jmod = abs(J*qw(n));
                    [c, gradc, ~, ~] = obj.localConcentrationInfo(R, dR, elm, e, obj.c0, obj.c0);
                    [~, epsilon] = obj.cf_asb.localDisplacementInfo(R, dR, elm, e, obj.d0);
                    g = topologicalInterp(c, obj.inter_geo);            
                    D0 = obj.cf_asb.D0;

                    num = num +Jmod*(D2*obj.fc(c) +dot(gradc,gradc));
                    den = den +Jmod*(D2/E0)*dot(g*D0*epsilon,epsilon);
                end
            end
            eta = num/den;
        end          
        %% Assembly
        function [K, M, F] = assembleElasticitySystem(obj)
            [gbi, elm, er, lm, qp, qw, n_quad, ndof, nel_dof, nel] = obj.cf_asb.preLoopParser;
            K = zeros(ndof);
            M = K;
            F = zeros(ndof,1);
            for e=1:nel
                k = zeros(nel_dof*obj.cf_asb.dimensions);
                m = zeros(nel_dof*obj.cf_asb.dimensions);
                f_e = zeros(nel_dof,obj.cf_asb.dimensions);
                for n=1:n_quad
                    q = qp(n,:);
                    [R, dR, J] = FastShape(obj.domain, q, gbi, elm, er, e);
                    Jmod = abs(J*qw(n));
                    [c, ~, ~, ~] = obj.localConcentrationInfo(R, dR, elm, e, obj.c0, obj.c0);
                    D0 = obj.cf_asb.D0;
                    g = topologicalInterp(c, obj.inter_geo);
                    B = elementStiffness(dR, obj.cf_asb.dimensions);
                    N = kron(R', eye(obj.cf_asb.dimensions));
                    
                    k = k+Jmod*B'*(g*D0)*B;
                    m = m+Jmod*obj.cf_asb.rho*g*(N'*N);
                    
                    x = obj.domain.evalPointFromQuadrature(q, er, e);
                    ff = obj.force_function(x);
                    for i=1:obj.cf_asb.dimensions
                        f_e(:,i) = f_e(:,i) +Jmod*R*ff(i);
                    end
                end
                idx = lm(:,e)';
                K(idx,idx) = K(idx,idx) +k;
                M(idx, idx) = M(idx,idx) +m;
                F(idx) = F(idx) +f_e(:);
            end
            b = obj.domain.extract_boundaries;
            normals = [-1 0 0; 1 0 0; 0 -1 0; 0 1 0];
            for i=1:numel(obj.sides)
                side = obj.sides(i);
                nn = normals(side,:);
                t = obj.trac{i};
                elements = b{side,2};
                [qp, qw] = obj.boundary_quad_rule(side);
                n_quad = length(qw);
                for j=1:numel(elements)
                    e = elements(j);
                    t_e = zeros(nel_dof, 3);
                    for n=1:n_quad
                        q = qp(n,:);
                        [R, dR, J] = FastShape(obj.domain, q, gbi, elm, er, e);
                        Jmod = abs(J*qw(n));
                        t_e = t_e +Jmod*dR.*t;
                    end
                    idx = lm(:,e)';
                    t_e = t_e(:,1:obj.cf_asb.dimensions);
                    F(idx) = F(idx) +t_e(:);
                end
            end
        end
        function d = solveDisplacement(obj)
            free_dofs = setdiff(1:length(obj.c0)*2,obj.clamped_dofs);
            [K, ~, F] = obj.assembleElasticitySystem;
            d = zeros(length(obj.c0)*2,1);
            d(free_dofs) = K(free_dofs,free_dofs)\F(free_dofs);
        end
        
        function residual = assembleResidual(obj, theta, option)
            [gbi, elm, er, lm, qp, qw, n_quad, ndof, nel_dof, nel] = obj.preLoopParser;
            residual = zeros(2*ndof,1);
            Rc = zeros(ndof,1);
            Rm = zeros(ndof,1);
            D1 = obj.D(1);
            D2 = obj.D(2);
            D3 = obj.D(3);
            for e=1:nel
                r_c = zeros(nel_dof,1);
                r_m = r_c;
                for n=1:n_quad
                    q = qp(n,:);
                    [R, dR, J] = FastShape(obj.domain, q, gbi, elm, er, e);
                    Jmod = abs(J*qw(n));
                    [c0, gradc0, mu0, gradmu0] = obj.localConcentrationInfo(R, dR, elm, e, obj.c0, obj.mu0);
                    [c1, gradc1, mu1, gradmu1] = obj.localConcentrationInfo(R, dR, elm, e, obj.c1, obj.mu1);
                    
                    gradmu = theta*gradmu0 +(1-theta)*gradmu1;
                    
                    r_c = r_c +Jmod*((c1 -c0)*R +obj.dt*(dR*gradmu'));
                    
                    [~, epsilon] = obj.cf_asb.localDisplacementInfo(R, dR, elm, e, obj.d0);
                    f0 = obj.dfdc(c0);
                    f1 = obj.dfdc(c1);
                    [~, dg0] = topologicalInterp(c0, obj.inter_geo);
                    [~, dg1] = topologicalInterp(c1, obj.inter_geo);
                    D0 = obj.cf_asb.D0/obj.cf_asb.young_modulus;
                    ep0 = dot(dg0*D0*epsilon,epsilon);
                    ep1 = dot(dg1*D0*epsilon,epsilon);
                    
                    mu = theta*mu0 +(1-theta)*mu1;
                    ep = theta*ep0 +(1-theta)*ep1;
                    df = theta*f0 +(1-theta)*f1;
                    gradc = theta*gradc0 +(1-theta)*gradc1;
                    
                    r_m = r_m+Jmod*((mu -D2*df +D3*ep)*R -D1*(dR*gradc'));
                end
                idx = lm(:,e)';
                residual(idx) = residual(idx) +r_c;
                residual(idx+ndof) = residual(idx+ndof) +r_m;
            end
        end
        
        function tangent = assembleTangent(obj, option)
            [gbi, elm, er, lm, qp, qw, n_quad, ndof, nel_dof, nel] = obj.preLoopParser;
            Kc = zeros(ndof);
            Km = Kc;
            Kcm = Kc;
            Kmc = Kc;
            D1 = obj.D(1);
            D2 = obj.D(2);
            D3 = obj.D(3);
            for e=1:nel
                k_cc = zeros(nel_dof);
                k_cm = k_cc;
                k_mc = k_cc;
                k_mm = k_cc;
                for n=1:n_quad
                    q = qp(n,:);
                    [R, dR, J] = FastShape(obj.domain, q, gbi, elm, er, e);
                    Jmod = abs(J*qw(n));
                    
                    [c, ~, ~, ~] = obj.localConcentrationInfo(R, dR, elm, e, obj.c0, obj.c0);
                    [~, epsilon] = obj.cf_asb.localDisplacementInfo(R, dR, elm, e, obj.d0);
                    f = obj.d2fdc(c);
                    [~, dg] = topologicalInterp(c, obj.inter_geo);
                    D0 = obj.cf_asb.D0/obj.cf_asb.young_modulus;
                    ep = dot(dg*D0*epsilon,epsilon);
                    
                    k_cc = k_cc +Jmod*(R*R');
                    k_cm = k_cm +Jmod*obj.dt*obj.theta*(dR*dR');
                    k_mc = k_mc +Jmod*obj.theta*((-D2*f +D3*ep)*(R*R') +D1*(dR*dR'));
                    k_mm = k_mm +Jmod*obj.theta*(R*R');
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
        function timeLoop(obj, option)
            obj.d0 = obj.solveDisplacement;
            t = 0;
            steps = 1;
            while (t < obj.t_max) && (steps < obj.max_timesteps)
                t = t+obj.dt;
                 [converged, n_steps, res_norm] = obj.timeStep(obj.newton_max_steps, option);
                if converged
                    obj.c0 = obj.c1;
                    obj.mu0 = obj.mu1;
                    obj.solution_array(:,steps) = obj.c1;
                    [~, eB, eI] = obj.computeEnergies;
                    obj.d0 = obj.solveDisplacement;
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
                    obj.mu1 = obj.mu0;
                    t = t-obj.dt;
                    obj.dt = 0.5*obj.dt;
                end
            end
            obj.timetable = obj.timetable(1:steps-1,:);
            obj.solution_array = obj.solution_array(:,1:steps-1);
        end
        
        function [converged, n_steps, res_norm] = timeStep(obj, max_corrections, option)
            ndof = length(obj.c0);
            converged = 0;
            % Predictor Step
            obj.c1 = obj.c0;
            obj.mu1 = obj.mu0;
            % Multicorrector Step
            for i=1:max_corrections+1
                residual = obj.assembleResidual(obj.theta);
                res_norm = norm(residual);
                if res_norm < obj.res_tol
                    converged = 1;
                    break
                end
                tangent = obj.assembleTangent(obj.theta);
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