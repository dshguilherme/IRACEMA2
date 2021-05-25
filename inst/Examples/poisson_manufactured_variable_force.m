%% Example 2: Poisson with method of manufactured solutions different function

% Solution: u = 3x^3 -0.5y^4 -xy
% -lap(u) = f -> f = 6y^2 -18x
% with Dirichlet Conditions = 3x^3 -0.5y^4 -xy

% Geometry
clearvars
close all
clc
a = 10;
error_norm = zeros(10,1);
refinements = ceil(logspace(0,2,10) +2);
for idx =1:10
    clearvars -except idx error_norm refinements a
    P1 = [0 0 0 1];
    P2 = [1 0 0 1];
    U = [0 0 1 1];
    line1 = Geometry(1, {U}, {P1, P2}, [1]);
    line2 = bs_translation(line1,[0 1 0]);
    domain = bs_ruled_surface(line1, line2);

    % Refinement
    Xi = linspace(0,1,refinements(idx));
    Xi = Xi(2:end-1);
    domain.degree_elevate(1,1);
    domain.degree_elevate(1,2);
    domain.knot_refine(Xi,1);
    domain.knot_refine(Xi,2);
    
    % Assembly
    asb = Poisson(1,"gauss",1,domain);
    K = asb.build_stiffness;
    f = @(x) 6*x(2)^2 -18*x(1);
    F = asb.variable_force(f);

    % Boundary Conditions
    boundaries = domain.extract_boundaries;
    cpoints = boundaries(:,2);
    cpoints = cell2mat(cpoints(:));
    cpoints = unique(cpoints);
    P = domain.points;
    dirichlet_cp = P(cpoints);
    u = @(x) 3*x(1)^3 -0.5*x(2)^4 -x(1)*x(2);
    g = cellfun(u,dirichlet_cp);
    [d, F, solution] = asb.dirichlet_linear_solve(K,F,g,cpoints);
    % h = solution.plot_solution(1);
    solution.asb = asb;
    error_norm(idx) = solution.l2_error_norm(u,1)
end

    loglog(1./refinements,error_norm,'LineWidth',2)
    title('Method of manufactured solutions : Poisson equation convergence','FontWeight','bold')
    xlabel('h size', 'FontWeight','bold')
    ylabel('L_2 error norm','FontWeight','bold')
    grid on
    set(gca,'FontSize',20)