classdef DedeMixed < CahnHilliardMixed & handle
   
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
         function obj = DedeMixed(coupled_field_assembler, trac, h0, ...
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
            obj.c1 = obj.c0;
            obj.d0 = zeros(length(obj.c0)*2,1);
            obj.d1 = obj.d0;
            free_dofs = setdiff(1:length(obj.d0),obj.clamped_dofs);
            [res, K] = obj.displacementResidual(1);
            while norm(res(free_dofs)) > 1e-4
                dd = K(free_dofs, free_dofs)\-res(free_dofs);
                obj.d0(free_dofs) = obj.d0(free_dofs)+dd;
                d0 = obj.d0;
                [res, K] = obj.displacementResidual(1);
            end
        end
        
        function mu0 = initial_potential(obj)
            [gbi, elm, er, lm, qp, qw, n_quad, ndof, nel_dof, nel] = obj.preLoopParser;
            F = zeros(ndof,1);
            K = zeros(ndof);
            for e=1:nel
                k = zeros(nel_dof);
                f = zeros(nel_dof,1);
                for n=1:n_quad
                    q = qp(n,:);
                    [R, dR, J] = FastShape(obj.domain, q, gbi, elm, er, e);
                    Jmod = abs(J*qw(n));
                    [c, gradc, ~, ~] = obj.localConcentrationInfo(R, dR, elm, e, obj.c0, obj.c0);
                    df = obj.dfdc(c);
                    [d, epsilon] = obj.cf_asb.localDisplacementInfo(R, dR, elm, e, obj.d0);
                    D0 = obj.cf_asb.D0;
                    g = topologicalInterp(c, obj.inter_geo);
                    k = k+Jmod*(R*R');
                    f = f +Jmod*(df -1.5*obj.eta*dot(g*D0*epsilon,epsilon)*R +...
                        obj.lambda*dR*gradc');
                end
                idx = lm(:,e)';
                F(idx) = F(idx) +f;
                K(idx, idx) = K(idx, idx) +k;
            end
            b = obj.domain.extract_boundaries;
            dofs = unique(cell2mat(b(:,4)));
            free_dofs = setdiff(1:ndof, dofs);
            mu0 = zeros(ndof,1);
            mu0(free_dofs) = K(free_dofs,free_dofs)\F(free_dofs);
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
        function residual = concentrationResidual(obj, theta)
            [gbi, elm, er, lm, qp, qw, n_quad, ndof, nel_dof, nel] = obj.preLoopParser;
            residual = zeros(ndof,1);
            for e=1:nel
                r_e = zeros(nel_dof,1);
                for n=1:n_quad
                    q = qp(n,:);
                    [R, dR, J] = FastShape(obj.domain, q, gbi, elm, er, e);
                    Jmod = abs(J*qw(n));
                    [c0, ~, ~, gradmu0] = obj.localConcentrationInfo(R, dR, elm, e, obj.c0, obj.mu0);
                    [c1, ~, ~, gradmu1] = obj.localConcentrationInfo(R, dR, elm, e, obj.c1, obj.mu1);
                    gradmu = theta*gradmu0 +(1-theta)*gradmu1;
                    M = 1;
                    
                    r_e = r_e +Jmod*((c1 -c0)*R +M*obj.dt*(dR*gradmu'));
                end
                idx = lm(:,e)';
                residual(idx) = residual(idx) +r_e;
            end
            residual = sparse(residual);
        end
        
        function residual = potentialResidual(obj, theta)
            [gbi, elm, er, lm, qp, qw, n_quad, ndof, nel_dof, nel] = obj.preLoopParser;
            residual = zeros(ndof,1);
            D2 = obj.D(2);
            D3 = obj.D(3);
            D1 = obj.D(1);
            for e=1:nel
                r_e = zeros(nel_dof,1);
                r_e0 = r_e;
                r_e1 = r_e;
                for n=1:n_quad
                    q = qp(n,:);
                    [R, dR, J] = FastShape(obj.domain, q, gbi, elm, er, e);
                    Jmod = abs(J*qw(n));
                    [c0, gradc0, mu0, ~] = obj.localConcentrationInfo(R, dR, elm, e, obj.c0, obj.mu0);
                    [c1, gradc1, mu1, ~] = obj.localConcentrationInfo(R, dR, elm, e, obj.c1, obj.mu1);

                    [~, epsilon0] = obj.cf_asb.localDisplacementInfo(R, dR, elm, e, obj.d0);
                    [~, epsilon1] = obj.cf_asb.localDisplacementInfo(R, dR, elm, e, obj.d1);

                    f0 = obj.dfdc(c0);
                    f1 = obj.dfdc(c1);
                    
                    [~, dg0] = topologicalInterp(c0, obj.inter_geo);
                    [~, dg1] = topologicalInterp(c1, obj.inter_geo);
                    D0 = obj.cf_asb.D0;
                    ep0 = dot(dg0*D0*epsilon0,epsilon0);
                    ep1 = dot(dg1*D0*epsilon1,epsilon1);
                    
                    r_e0 = r_e0 +Jmod*((mu0 -f0 +obj.eta*ep0)*R -obj.lambda*(dR*gradc0'));
                    r_e1 = r_e1 +Jmod*((mu1 -D2*f1 +D3*ep1)*R -D1*(dR*gradc1'));
                    r_e = r_e+ theta*r_e0 +(1-theta)*r_e1;
                end
                idx = lm(:,e)';
                residual(idx) = residual(idx) +r_e1;
            end
        end
        
        function [residual, K] = displacementResidual(obj, theta)
            [gbi, elm, er, lm, qp, qw, n_quad, ndof, nel_dof, nel] = obj.cf_asb.preLoopParser;
            residual = zeros(ndof,1);
            K = zeros(ndof);
            K1 = zeros(ndof);
            F = zeros(ndof,1);
            d = obj.cf_asb.dimensions;
            D4 = obj.D(4);
            D5 = obj.D(5);
            for e=1:nel
                k = zeros(nel_dof*d);
                k0 = k;
                k1 = k;
                f_e = zeros(nel_dof,d);
                for n=1:n_quad
                    q = qp(n,:);
                    [R, dR, J] = FastShape(obj.domain, q, gbi, elm, er, e);
                    Jmod = abs(J*qw(n));
                    [c0, ~, ~, ~] = obj.localConcentrationInfo(R, dR, elm, e, obj.c0, obj.c0);
                    [c1, ~, ~, ~] = obj.localConcentrationInfo(R, dR, elm, e, obj.c1, obj.c1);

                    [g0, ~] = topologicalInterp(c0, obj.inter_geo);
                    [g1, ~] = topologicalInterp(c1, obj.inter_geo);
                    B = elementStiffness(dR, d);
                    D0 = obj.cf_asb.D0;
                    k0 = k0+Jmod*2*(B'*g0*D0*B);
                    k1 = k1+Jmod*2*(B'*g1*D0*B);
                    x = obj.domain.evalPointFromQuadrature(q, er, e);
                    ff = obj.force_function(x);
                    for i=1:d
                        f_e(:,i) = f_e(:,i) +Jmod*R*ff(i);
                    end
                end
                idx = lm(:,e)';
                F(idx) = F(idx) +f_e(:);
                K(idx,idx) = K(idx,idx) +k0;
                K1(idx, idx) = K1(idx, idx) +k1;
            end
            K = sparse(K);
            K1 = sparse(K1);
            F = sparse(F);
            residual = (K1*obj.d1 -F);
        end
        
        function tangent = assembleTangents(obj, K)
            [gbi, elm, er, lm, qp, qw, n_quad, ndof, nel_dof, nel] = obj.preLoopParser;
            Kcc = zeros(ndof);
            Kcm = Kcc;
            Kmc = Kcc;
            Kmm = Kcc;
            Kmd = zeros(ndof, obj.cf_asb.dimensions*ndof);
            Kcd = Kmd;
            Kdc = Kmd';
            Kdm = Kdc;
            D1 = obj.D(1);
            D2 = obj.D(2);
            D3 = obj.D(3);
            for e=1:nel
                kcc = zeros(nel_dof);
                kcm = kcc;
                kmc = kcc;
                kmm = kcc;
                kmd = zeros(nel_dof,nel_dof*obj.cf_asb.dimensions);
                kdc = kmd';
                for n=1:n_quad
                    q = qp(n,:);
                    [R, dR, J] = FastShape(obj.domain, q, gbi, elm, er, e);
                    Jmod = abs(J*qw(n));
                    [c, ~, ~, ~] = obj.localConcentrationInfo(R, dR, elm, e, obj.c0, obj.mu0);
                    [~, epsilon] = obj.cf_asb.localDisplacementInfo(R, dR, elm, e, obj.d0);
                    M = obj.mobility;
                    f = obj.d2fdc(c);
                    D0 = obj.cf_asb.D0/obj.cf_asb.young_modulus;
                    [~, dg] = topologicalInterp(c, obj.inter_geo);
                    B = elementStiffness(dR, obj.cf_asb.dimensions);
                    N = B'*dg*D0*epsilon;
                    dep = dot(dg*D0*epsilon,epsilon);
                    
                    kcc = kcc+Jmod*(R*R');
                    kcm = kcm+Jmod*M*obj.dt*obj.theta*(dR*dR');
                    kmc = kmc -Jmod*((D2*f +D3*dep*(R*R')) +D1*(dR*dR'));
                    kmm = kmm +Jmod*(R*R');
                    kmd = kmd +Jmod*D3*(R*N');
                    kdc = kdc +Jmod*(N*R');
                end
                idx = lm(:,e)';
                Kcc(idx,idx) = Kcc(idx,idx)+ kcc;
                Kcm(idx,idx) = Kcm(idx,idx)+ kcm;
                Kmm(idx,idx) = Kmm(idx,idx)+ kmm;
                Kmc(idx,idx)= Kmc(idx,idx)+ kmc;
                Kmd(idx, [idx ndof+idx]) = Kmd(idx, [idx ndof+idx]) +kmd;
                Kdc([idx; ndof+idx], idx) = Kdc([idx; ndof+idx], idx) +kdc;
            end

            tangent = sparse([Kcc, Kcm, Kcd;
                       Kmc, Kmm, Kmd;
                       Kdc, Kdm, K]);
                    
        end
        %% Time Integration Functions
        function [converged, n_steps, res_norm, Rc_BE, Rc_CN] = timeStep(obj, max_corrections)
            ndof = length(obj.c0);
            converged = 0;
            % Predictor Step
            obj.c1 = obj.c0;
            obj.mu1 = obj.mu0;
            obj.d1 = obj.d0;
            % Boundary condition indexes
            skip = 2*length(obj.c0);
            fd = setdiff(1:length(obj.d0),obj.clamped_dofs);
            % Multicorrector Step
            for i=1:max_corrections+1
                Rc_CN = obj.concentrationResidual(obj.theta);
                Rm_CN = obj.potentialResidual(obj.theta);
                [Rd_CN, K] = obj.displacementResidual(obj.theta);
                residual = [Rc_CN;Rm_CN;Rd_CN];
                norm_c = norm(Rc_CN);
                norm_mu = norm(Rm_CN);
                norm_d = norm(Rd_CN(fd));
                res_norm = norm([norm_c, norm_mu, norm_d]);
                if res_norm < obj.res_tol
                    converged = 1;
                    break
                end
                tangent = obj.assembleTangents(K);
                u = zeros(size(residual));
                dofs = length(obj.c0)*2+obj.clamped_dofs;
                free_dofs = setdiff(1:length(residual), dofs);
                switch obj.mode
                    case "Newton"
                        u(dofs) = 0;
                        u(free_dofs) = tangent(free_dofs, free_dofs)\-(residual(free_dofs));
                    case "GMRES"
                        tangent(dofs,:) = [];
                        tangent(:,dofs) = [];
                        residual(dofs) = [];
%                         [L, U] = ilu(tangent,struct('type','ilutp','droptol',1e-6));
                        [uu, flag, ~, ~, ~] = gmres(tangent,-residual,[],obj.gmres_tol,obj.gmres_maxit);
                        if flag
                            converged = 0;
                            break
                        end
                        u(free_dofs) = uu;
                        u(dofs) = 0;
                end
                obj.c1 = obj.c1 +u(1:ndof);
                obj.mu1 = obj.mu1 +u(ndof+1:2*ndof);
                obj.d1 = obj.d1 +u(2*ndof+1:end);
            end
            n_steps = i;
            Rc_BE = obj.concentrationResidual(1);
        end
        
        function timeLoop(obj)
            t = 0;
            % Initialize vectors
            steps = 1;
            while (t < obj.t_max) && (steps < obj.max_timesteps)
                t = t+obj.dt;
                [converged, n_steps, res_norm, R_BE, R_CN] = obj.timeStep(obj.newton_max_steps);
                r = norm(R_CN - R_BE)/norm(R_CN);
                if converged && (n_steps < 8) && (r < 1.2*obj.res_tol)
                    obj.c0 = obj.c1;
                    obj.mu0 = obj.mu1;
                    obj.d0 = obj.d1;
                    obj.solution_array(:,steps) = obj.c1;
                    [~, eB, eI] = obj.computeEnergies;
                    [eP, eK] = obj.cf_asb.computeEnergies(obj.frequency, obj.d0);
                    eT = eB+eI+eP;
                    obj.timetable(steps,:) = [t, obj.dt, eT, eB, eI, eP, res_norm];
                    fprintf("Step: %i | Time: %.2e | dt: %.2e | Residual Norm: %.2e \n", steps, t, obj.dt , res_norm)
                    r = 1.2*norm(R_CN - R_BE)/norm(R_CN);
                    obj.dt = obj.dt*((obj.res_tol/r)^0.25);
                    steps = steps+1;
                else
                    obj.c1 = obj.c0;
                    obj.mu1 = obj.mu0;
                    t = t-obj.dt;
                    obj.dt = 0.5*obj.dt;
                end
            end
            obj.timetable = obj.timetable(1:steps-1,:);
        end
        
     
        
        %% Auxiliary Functions
        function plotEvolution(obj)
                t_array = obj.timetable(:,1);
                dt_array = obj.timetable(:,2);
                eT = obj.timetable(:,3);
                eB = obj.timetable(:,4);
                eI = obj.timetable(:,5);
                eP = obj.timetable(:,6);
%                 if option == "dynamic"
%                     eK = obj.timetable(:,7);
%                 end
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
%                 if option == "dynamic"
%                     plot(t_array,eK, 'LineWidth', 3)
%                     legend('Total Energy', 'Bulk Energy', 'Interface Energy', 'Elastic Energy', 'Kinetic Energy')
%                 else
                    legend('Total Energy', 'Bulk Energy', 'Interface Energy', 'Elastic Energy')
%                 end
                set(gca, 'FontSize', 18)
                ylabel('Energy [J]', 'FontWeight', 'bold', 'FontSize', 23)
                xlabel('time [s]', 'FontWeight', 'bold', 'FontSize', 23)
                title('Free Energy Functional', 'FontWeight', 'bold', 'FontSize', 23)
                grid on
           end
        
        
        
    end

end