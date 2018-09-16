%% Ejemplo térmica
clear; close all; clc
L = 2000; %mm
E = [70000 100000 100000]; %MPa
alfa = [23 20 20]*10^-6; 
A = [12 6 6]*10^2; %mm^2
tita = -10;
FT = abs(tita*alfa(1)*E(1)*A(1));
Sig0=-tita*alfa(1)*E(1);
%% Nodos
cnod = [0 L 2*L 2*L];
nnod = length(cnod);
elem = [1 2
        2 3
        2 4];
nelem = size(elem,1);
%% Dof
dofpornodo = 1;
dof = [1
       2
       3
       4];
%% Matriz global
Kglobal = zeros(4);
for e = 1:nelem
    vec = cnod(elem(e,2))-cnod(elem(e,1));
    long = norm(vec);
    Klocal = A(e)*E(e)/long*[1 -1; -1 1];
    dofr = [dof(elem(e,1)) dof(elem(e,2))];
    Kglobal(dofr,dofr) = Kglobal(dofr,dofr)+Klocal;
end
%% Cargas térmicas
P = [FT -FT 0 0];
Pred = P(2);
Kglobalred = Kglobal(2,2);
%% Solver
Dred = Pred/Kglobalred;
D = [0 Dred 0 0]';
P = Kglobal*D-P';
%% Tensiones
Sigma=[-P(1)/A(1) P(3)/A(2) P(4)/A(3)]'