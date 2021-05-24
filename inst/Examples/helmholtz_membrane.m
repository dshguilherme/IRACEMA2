clearvars
clc
close all
%% Example 6: Harmonic Analysis of a Membrane
L = pi;
h = exp(1);
P1 = [0 0 0 1];
P2 = [L 0 0 1];

line = Geometry(1,{[0 0 1 1]}, {P1,P2},1);
line2 = bs_translation(line,[0 h 0]);

domain = bs_ruled_surface(line,line2);

rho = 1e3;
alpha = 1e3;

% Refinement
domain.degree_elevate(1,1);
domain.degree_elevate(1,2);

Xi = linspace(0,1,21);
Xi = Xi(2:end-1);
domain.knot_refine(Xi,1);
domain.knot_refine(Xi,2);

% Assembly
asb = Helmholtz(rho,alpha,"gauss",1,domain);
[K, M] = asb.build_matrices;

% Boundary Conditions
clamp_dofs = domain.extract_boundaries;
clamp_dofs = clamp_dofs(:,2);
clamp_dofs = unique(cell2mat(clamp_dofs(:)));

% Solution
[vecs, omega, solution] = asb.eigensolve(K,M,10,clamp_dofs);

z = 1:numel(omega);
y = zeros(numel(omega),1);

for i=1:numel(Xi)+3
    for j=1:numel(Xi)+3
        wn(i,j) = pi*sqrt((i/L)^2 +(j/h)^2);
    end
end
wn = sort(wn(:));
y = omega./(wn(1:numel(y)));
N = (1:numel(y))/numel(y);
plot(N,y);
grid on
set(gca,'FontSize',18)
ylabel('\omega^h / \omega_n', 'FontWeight','bold','FontSize',20)
xlabel('n/N', 'FontWeight','bold','FontSize',20)


