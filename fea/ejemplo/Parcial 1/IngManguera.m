%% Ej de desp conocidos con P conocidos
clear; close all; clc
%% Datos
E = 12500;
h = 50;
b = 300;
I = b*h^3/12;
L = 1000;
M2 = 1000000; %N.mm
U3 = -10; %mm
%% Matriz
Kglobal = zeros(6);
Y4=2*E*I/L;
Y3=Y4*2;
Y2=Y4*3/L;
Y1=Y2*2/L;
Kviga=[Y1 Y2 -Y1 Y2
       Y2 Y3 -Y2 Y4
      -Y1 -Y2 Y1 -Y2
       Y2 Y4 -Y2 Y3];
Kglobal(1:4,1:4) = Kglobal(1:4,1:4)+Kviga;
Kglobal(3:6,3:6) = Kglobal(3:6,3:6)+Kviga;
%% BC
% El desplazamiento conocido es uan condicion de borde
fijo = 1&[1 1 1 0 1 0];
libre = ~fijo;
%% Cargas conocidas
P = [0 0 0 M2 0 0]';
%% Desplazamientos conocidos
D = [0 0 0 0 U3 0]';
%% Reduccion
Kglobalx = Kglobal(libre,libre);
Kglobalcx = Kglobal(libre,fijo);
Kglobalc = Kglobal(fijo,fijo);
Pc = P(libre);
Dc = D(fijo);
%% Solver
Dx = Kglobalx\(Pc-Kglobalcx*Dc);
D(libre) = Dx;
P = Kglobal*D;
disp(strcat(['La poronga del Ing. Maguera pesa: ',num2str(-P(5)/10-80),' Kg, flojo pedazo.']))