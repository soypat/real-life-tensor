%% Ej 1
%clear; close all; clc
% Por simetría, el nodo en el medio no va a desplazarse por su condicion de
% brode ni rotarse por la simetría. Por esto, invoco simetría y supongo en
% el centro una condición de borde de empotramiento.
function [Dpunta]=DespPunta(alfa)
%% Datos
a = 15; %mm
I = a^4/12; %mm^4
E = 210000; %MPa
L = 1000; %mm
k = 100; %N/mm
q = -1; %N/mm
%alfa = 0.2599;
L1 = alfa*L/2;
L2 = (1-alfa)*L/2;
%% Nodos y elementos
cnod = [0 0
        L1 0
        L1+L2 0];
nnod = size(cnod,1);
elem = [1 2
        2 3];
nelem = size(elem,1);
%% Dof
dofpornodo = 2;
dof = reshape((1:1:dofpornodo*nnod)',dofpornodo,nnod)';
%% Matriz de rigidez
Kglobal = zeros(dofpornodo*nnod);
memlong = zeros(1,nelem);
for e = 1:nelem
    vec = cnod(elem(e,2),:) - cnod(elem(e,1),:);
    long = norm(vec);
    memlong(e) = long;
    %vd = vec/long;
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
Kglobal(1,1) = Kglobal(1,1)+k; %Inclusion del resorte
%% BC
fijo = [0 0 1 0 1 1]';
libre = ~fijo;
%% Cargas
P = [-q*memlong(1)/2 -q*memlong(1)^2/12 
    -q*memlong(1)/2-q*memlong(2)/2 q*memlong(1)^2/12-q*memlong(2)^2/12
    -q*memlong(2)/2 q*memlong(2)^2/12];
P = reshape(P',[],1);
%% Reductor
Kglobalred = Kglobal(libre,libre);
Pred = P(libre);
%% Solver
Dred = Kglobalred\Pred;
D = [Dred(1:2)' 0 Dred(3) 0 0]';
Dpunta = D(1);
end