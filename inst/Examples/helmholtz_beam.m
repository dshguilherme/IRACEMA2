clearvars
clc
%% Example 6: Helmholtz Equation of a Beam
rho = 1e3;
alpha = 1e3;
% Geometry
L = 1.0;
P1 = [0 0 0];
P2 = [L 0 0];

domain = bs_line(P1,P2);
domain2 = bs_line(P1,P2);

% Refinement
% IGA:
domain.degree_elevate(1,1);
Xi = linspace(0,1,1000);
Xi = Xi(2:end-1);
domain.knot_refine(Xi,1);

% Assembly
asb = Helmholtz(rho,alpha,"gauss",1,domain);
[K, M] = asb.build_matrices;
% BCs
[K_c, M_c] = asb.clamp_boundary(K,M);
% Solution
[~, omega] = eigs(K_c,M_c,numel(Xi)+1,'sm');
omega = sqrt(diag(omega));

n = 1:numel(omega);
y = zeros(numel(omega),1);
N = y;
for i=n
    wn = n(i)*pi;
    y(i) = omega(i)./wn;
    N(i) = n(i)/length(n);
end
figure(1)
plot(N,y);
ylim([1 1.15])

