%% Ej 1
clear; close all; clc
%% Alfa
n = 100;
alfa = linspace(0.1,0.9,n);
D = zeros(n,1);
for i = 1:n
    D(i) = DespPunta(alfa(i));
end
minimo = min(abs(D));
m=find(D==minimo);
alfa = linspace(alfa(m-1),alfa(m+1),n);
D = zeros(n,1);
for i = 1:n
    D(i) = DespPunta(alfa(i));
end
minimo = min(abs(D));
m=find(abs(D)==minimo);
disp(alfa(m))