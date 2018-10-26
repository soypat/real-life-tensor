%% Ej 1 solución exacta
clear; close all; clc
%% Datos
c = 10;
L = 100;
w = 1;
P = 80;
E = 1000;
nu = .25;

%% Solución exacta
syms x y real
I = (2*c)^3*w/12;
G = E/2/(1+nu);
u = P/(2*I)*(-(x^2-L^2)/E*y-nu*y*(y^2-c^2)/(3*E)+y*(y^2-c^2)/(3*G));
v = P/I*(nu*x*y^2/(2*E)+(x^3-L^3)/(6*E)-(L^2/(2*E)+nu*c^2/(6*E)+c^2/(3*G))*(x-L));
e11 = diff(u,x);
e22 = diff(v,y);
e12 = (diff(u,y)+diff(v,x))/2;
e11 = matlabFunction(e11);
e22 = matlabFunction(e22);
e12 = matlabFunction(e12);
[X,Y] = meshgrid(linspace(0,L,101),linspace(-c,c,21));
epsilon = zeros(21,101,3);
epsilon(:,:,1) = e11(X,Y);
epsilon(:,:,2) = e22(X,Y);
epsilon(:,:,3) = e12(Y);

