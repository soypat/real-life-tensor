%% Patch test Q4
clear; close all; clc
%% Datos
E = 210E9;
nu = 0.3;
eleT = 'Q4';
t = 1;
f = 1000;
%% Nodos
nod = [0 0
       1 0
       2 0
       0 1
       1.2 0.7
       2 1
       0 2
       1 2
       2 2];
nnod = size(nod,1);
elem = [1 2 5 4
        2 3 6 5
        4 5 8 7
        5 6 9 8];
nelem = size(elem,1);

%% dof
dofpornodo = 2;
nodporelem = 4;
doftot = dofpornodo*nnod;
dof = reshape((1:doftot)',dofpornodo,nnod)';

%% Constitutiva
% Plane stress
C = [1 nu 0
     nu 1 0
     0 0 1-nu]*E/(1-nu^2);

%% Matriz
[wpg, upg, npg] = gauss([2 2]); 
Kglobal = zeros(doftot);
for e = 1:nelem
    nodelem = nod(elem(e,:),:);
    Ke = zeros(nodporelem*dofpornodo);
    for ipg = 1:npg
        ksi = upg(ipg,1);
        eta = upg(ipg,2);
        dNke = shapefunsder([ksi eta],eleT);
        J = dNke*nodelem;
        dNxy = J\dNke;
        B = zeros(size(C,2),nodporelem*dofpornodo);
        B(1,1:2:7) = dNxy(1,:);
        B(2,2:2:8) = dNxy(2,:);
        B(3,1:2:7) = dNxy(2,:);
        B(3,2:2:8) = dNxy(1,:);
        Ke = Ke + B'*C*B*t*wpg(ipg)*det(J);
    end
    dofs = reshape(dof(elem(e,:),:)',[],1);
    Kglobal(dofs,dofs) = Kglobal(dofs,dofs) + Ke;
end

%% Cargas
P = zeros(doftot/2,2);
P(nod(:,1)==2,1) = f*[1 2 1]';
P = reshape(P',[],1);

%% BC
fijo = zeros(doftot/2,2);
fijo(1,:) = true;
fijo([4 7],1) = true;
fijo = reshape(fijo',[],1);
libre = ~fijo;

%% Solver
Dred = Kglobal(libre,libre)\P(libre);
D = zeros(doftot,1);
D(libre) = Dred;

%% Deformación
strain = zeros(nelem*nodporelem,3);
[wpg, upg, npg] = gauss([2 2]);
rsextra = [-1 -1
           -1 1
           1 -1
           1 1]*sqrt(3);
for e = 1:nelem
    nodelem = nod(elem(e,:),:);
    dofs = reshape(dof(elem(e,:),:)',[],1);
    for ipg = 1:npg
        pgstrain = zeros(4,3)
        ksi = upg(ipg,1);
        eta = upg(ipg,2);
        dNke = shapefunsder([ksi eta],eleT);
        J = dNke*nodelem;
        dNxy = J\dNke;
        B = zeros(size(C,2),nodporelem*dofpornodo);
        B(1,1:2:7) = dNxy(1,:);
        B(2,2:2:8) = dNxy(2,:);
        B(3,1:2:7) = dNxy(2,:);
        B(3,2:2:8) = dNxy(1,:);
        pgstrain(ipg,:) = B*D(dofs);
        
        
    end
end






















