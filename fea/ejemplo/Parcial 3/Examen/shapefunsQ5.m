%% Funciones de forma para Q5
clear; close all; clc

syms x y real
X = [1   x   y   x*y   x^2*y^2];
A = [1  -1  -1    1       1
     1   1  -1   -1       1
     1   1   1    1       1
     1  -1   1   -1       1
     1   0   0    0       0];
N = X/A;
dNx = diff(N,x);
dNy = diff(N,y);
dN = [dNx; dNy];
% sum(N) % == 1
% sum(sum(dN)) % == 0