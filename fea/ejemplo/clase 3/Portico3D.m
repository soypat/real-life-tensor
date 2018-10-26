%% Pórtico 3D
clear; close all; clc
%% Datos
E = 210000; %Mpa
nu = 0.3;
G = E/2/(1+nu);
Fx = -5000; %N
Fz = 2000; %N
q = 10; %N/mm
tita = 20; %'C de la columna 1
alfa = 1.2E-5; 
% Columnas IPE 300
Ac = 53.8E2; %mm^2 Area perfil columnas
Izc = 8360E4; %mm^4 Momento de iniercia en y local
Iyc = 604E4; %mm^4 Momento de inercia en z local
Jxc = 15.6E4; %mm^4 Módulo de torsión
% Vigas IPE 160
Av = 20.1E2; %mm^2 Area perfil vigas
Izv = 869E4; %mm^4 Momento de iniercia en y local
Iyv = 68.3E4; %mm^4 Momento de inercia en z local
Jxv = 2.82E4; %mm^4 Módulo de torsión
%% Nodos y elementos
cnod = [0 0 0
        0 3000 0
        2000 4000 0
        4000 3000 0
        4000 0 0];
nnod = size(cnod,1);
aux = [-1 0 0];
elem = [1 2
        2 3
        3 4
        4 5];
nelem = size(elem,1);
A = [Ac Av Av Ac];
Iz = [Izc Izv Izv Izc];
Iy = [Iyc Iyv Iyv Iyc];
J = [Jxc Jxv Jxv Jxc];
%% Matriz global
dofpornodo = 6; %Dof por nodo
dof = reshape((1:1:dofpornodo*nnod)',dofpornodo,nnod)';
Kglobal = zeros(dofpornodo*nnod);
memlong = zeros(nelem,1);
memvecdirec = zeros(nelem,3);
memT = zeros(12,12,nelem);
for e = 1:nelem
    vec = cnod(elem(e,2),:)-cnod(elem(e,1),:);
    long = norm(vec);
    vecdirec = vec/long;
    memlong(e) = long;
    memvecdirec(e,:) = vecdirec;
%     Lambda = [vecdirec(1:2) 0; -vecdirec(2) vecdirec(1) 0; 0 0 1];
%     T = blkdiag(Lambda,Lambda,Lambda,Lambda);
%     memT(:,:,e) = T;
    X = A(e)*E/long;
    Y4 = 2*E*Iz(e)/long;
    Y3 = Y4*2;
    Y2 = 3*Y4/long;
    Y1 = 2*Y2/long;
    Z4 = 2*E*Iy(e)/long;
    Z3 = Z4*2;
    Z2 = 3*Z4/long;
    Z1 = 2*Z2/long;
    S = G*J(e)/long;
    Klocal = [X 0 0 0 0 0 -X 0 0 0 0 0
              0 Y1 0 0 0 Y2 0 -Y1 0 0 0 Y2
              0 0 Z1 0 -Z2 0 0 0 -Z1 0 -Z2 0
              0 0 0 S 0 0 0 0 0 -S 0 0
              0 0 -Z2 0 Z3 0 0 0 Z2 0 Z4 0
              0 Y2 0 0 0 Y3 0 -Y2 0 0 0 Y4
              -X 0 0 0 0 0 X 0 0 0 0 0
              0 -Y1 0 0 0 -Y2 0 Y1 0 0 0 -Y2
              0 0 -Z1 0 Z2 0 0 0 Z1 0 Z2 0
              0 0 0 -S 0 0 0 0 0 S 0 0
              0 0 -Z2 0 Z4 0 0 0 Z2 0 Z3 0
              0 Y2 0 0 0 Y4 0 -Y2 0 0 0 Y3];
    Kelem = rotar(Klocal,cnod(elem(e,1),:),cnod(elem(e,2),:),aux);
    dofr = [dof(elem(e,1),:) dof(elem(e,2),:)];
    Kglobal(dofr,dofr) = Kglobal(dofr,dofr)+Kelem;
end
%% BC
fijo = [ones(1,6) zeros(1,18) ones(1,6)];
libre = ~fijo;
%% Cargas
FT = alfa*tita*A(1)*E; % Gradientes de Temperatura
Pm = [zeros(1,6) 
      0 FT Fz 0 0 0
      -q*memlong(3)/2*-memvecdirec(3,2) -q*memlong(3)/2*memvecdirec(3,1) Fz 0 0 -q*memlong(3)^2/12
      -q*memlong(3)/2*-memvecdirec(3,2)+Fx -q*memlong(3)/2*memvecdirec(3,1) 0 0 0 q*memlong(3)^2/12
      zeros(1,6)];
P = reshape(Pm',[],1);
%% Reducción
Kglobalred = Kglobal(libre,libre);
Pred = P(libre);
%% Solver
Dred = Kglobalred\Pred;
D = [zeros(1,6) Dred' zeros(1,6)];
