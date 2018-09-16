%% Estructura de datos geométricos
clear all; close all; clc
L=60; A=2; E=30*10^6;
neletot=input('Ingrese cantidad de elementos: ');
des=sparse(1,neletot);
for nele=1:neletot
    des(nele)=linear_load(L,A,E,nele);
end
plot(1:nele,des);
grid on; xlabel('N de Elementos'); ylabel('Desplazamiento en x=L (in)')