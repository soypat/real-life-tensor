%% Ejercicio re crazy
clear ; close all; clc
%% Datos
a = 50;      %mm
I = a^4/12;  %mm^4
L = 1000;    %mm
E = 210000;  %MPa
k = 1000;    %N/mm
q = 10;      %N/mm
%% Nodos y elementos
cnod = [0 0
        L 0
        2*L 0
        0 0
        L 0];
nnod = 5;
elem = [1 2
        2 3
        4 5];
nelem = 3;
%% Dof
dofpornodo = 2;
dof = reshape((1:1:dofpornodo*nnod)',dofpornodo,nnod)';
%% Matriz de rigidez
Kglobal = zeros(dofpornodo*nnod);
long = L;
for e = 1:nelem
    Y4 = 2*E*I/long;
    Y3 = Y4*2;
    Y2 = Y4*3/long;
    Y1 = Y2*2/long;
    Kviga = [Y1 Y2 -Y1 Y2
             Y2 Y3 -Y2 Y4
            -Y1 -Y2 Y1 -Y2
             Y2 Y4 -Y2 Y3];
    dofr = [dof(elem(e,1),:) dof(elem(e,2),:)];
    Kglobal(dofr,dofr) = Kglobal(dofr,dofr)+Kviga;
end
% Resorte 2 con 5
Kglobal([3 9],[3 9]) = Kglobal([3 9],[3 9])+k*[1 -1;-1 1];
% Resorte 3
Kglobal(5,5) = Kglobal(5,5)+k;
%% BC
fijo = [1 1 0 0 0 0 1 1 0 0]';
libre = ~fijo;
%% Carga
P = [0 0 0 0 0 0 -q*L/2 -q*L^2/12 -q*L/2 q*L^2/12]';
%% Reducción
Kglobalred = Kglobal(libre,libre);
Pred = P(libre);
%% Solver
Dred = Kglobalred\Pred;
D = [0 0 Dred(1:4)' 0 0 Dred(5:6)']';
Dpornodo = reshape(D,dofpornodo,nnod)'