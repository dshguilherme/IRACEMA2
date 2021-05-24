%% Example 2: Poisson with method of manufactured solutions different function

% Solution: u = 2^(4*a)*(x^a)*((1-x)^a)*(y^a)*(1-y)^a
% -lap(u) = f -> f = 5*(pi^2)*sen(pi*x)*sen(pi*y)
% with Dirichlet Conditions = 0

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
    fx1 = @(x) (2^(4*a))*(x(2)^a)*((1-x(2))^a);
    fx2 = @(x) (x(1)^a)*a*((a-1)*(1-x(1))^(a-2)) -a*((1-x(1))^(a-1)*a*x(1)^(a-1)) ...
        +a*(x(1)^(a-1))*(-a*((1-x(1))^(a-1))) +((1-x(1))^a)*a*(a-1)*(x(1)^(a-2));
    fy1 = @(x) (2^(4*a))*(x(1)^a)*((1-x(1))^a);
    fy2 = @(x) (x(2)^a)*a*((a-1)*(1-x(2))^(a-2)) -a*((1-x(2))^(a-1)*a*x(2)^(a-1)) ...
        +a*(x(2)^(a-1))*(-a*((1-x(2))^(a-1))) +((1-x(2))^a)*a*(a-1)*(x(2)^(a-2));
    f = @(x) -(fx1(x)*fx2(x) +fy1(x)*fy2(x));
    F = asb.variable_force(f);

    % Boundary Conditions
    boundaries = domain.extract_boundaries;
    cpoints = boundaries(:,2);
    cpoints = cell2mat(cpoints(:));
    cpoints = unique(cpoints);
    P = domain.points;
    dirichlet_cp = P(cpoints);
    u = @(x) (2^(4*a))*(x(1)^a)*((1-x(1))^a)*(x(2)^a)*((1-x(2))^a);
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