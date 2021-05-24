close all
clearvars
clc

%% Geometric and Material data

P1 = [-1 0 0];
P2 = [1 0 0];

line1 = bs_line(P1,P2);
line2 = bs_translation(line1,[0 1 0]);

domain = bs_ruled_surface(line1,line2);

rho = 1;
kappa = 1;