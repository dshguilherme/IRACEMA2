clearvars
close all
clc
% Geometry
P1 = [0 0 0 1];
P2 = [1 0 0 1];
U = [0 0 1 1];
line1 = Geometry(1, {U}, {P1, P2}, [1]);
line2 = bs_translation(line1,[0 1 0]);
domain = bs_ruled_surface(line1, line2);

% Material properties
lambda = 0.0005;
theta = 3/2;
diffusivity = 50;

% Refinement
Xi = linspace(0,1,10);
Xi = Xi(2:end-1);
domain.degree_elevate(1,1);
domain.degree_elevate(1,2);
domain.knot_refine(Xi,1);
domain.knot_refine(Xi,2);

% Initial condition
initial_c = 0.5 +0.05.*randn(numel(domain.points),1);

% Assembler
max_steps = 1000;
asb = CahnHilliard(lambda, 1, domain, 5e-6, 1, 500);


%% Solving
t = 0;
T = 1;
count = 0;
c = asb.concentration(:,1);
while (t < T) && (count < max_steps)
    t = t+asb.dt;
    count = t/asb.dt;
    asb.concentration(:,count) = c;
    % Predictor Step
    c0 = asb.concentration(:,count);
    c = c0;
    mu = asb.chemicalPotential(c0);
    % Multicorrector Steps
    for i=1:100
        c = c +(c-c0);
        [K, R] = asb.assembleMixedSystem(c, c0);
        if norm(R) < 1e-4
            break
        end
        deltac = K\(-R);
        c = c+ asb.dt*deltac(1:length(c));
    end
%     [eT, eB, eI] = asb.computeEnergies(c);
%     obj.timetable(count,:) = [time, dt, eT, eB, eI];
    disp("Step:");
    disp(count);
%     disp("Total Energy:");
%     disp(eT);
end
        
  