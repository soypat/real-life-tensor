%% Práctica 8 el culo te abrocho
clear; close all; clc
%% Datos
E = 5E9; %[Pa]
nu = 0.3;
lam = E*nu/(1+nu)/(1-2*nu);
mu = E/2/(1+nu);
t = 1; % espesor

%% Nodos y elementos
nod = load('nodos.txt');
nod = nod(:,[2 3]);
nnod = size(nod,1);
elem = load('elementos.txt');
elem = elem(:,2:5);
nelem = size(elem,1);

%% DOF
dofpornodo = 2;
nodporelem = 4;
doftot = dofpornodo*nnod;
dof = reshape(1:doftot,dofpornodo,nnod)';

%% Matriz constitutiva
C = [1 nu 0
     nu 1 0
     0 0 1-nu]*E/(1-nu^2);
 
%% Puntos de gauss
a = 1/sqrt(3);
epg = [-a -a
       a -a
       a a
       -a a];
npg = size(epg,1);

%% Matriz
Kglobal = zeros(doftot);
for e = 1:nelem
    Ke = zeros(dofpornodo*nodporelem);
    nodelem = nod(elem(e,:),:);
    for ig = 1:npg
        xi = epg(ig,1);
        eta = epg(ig,2);
        dNxieta = [-(1-eta) 1-eta 1+eta -(1+eta)
                   -(1-xi) -(1+xi) 1+xi 1-xi]/4;
        J = dNxieta*nodelem;
        dNxy = J\dNxieta;
        B = zeros(size(C,2),dofpornodo*nodporelem);
        B(1,1:2:7) = dNxy(1,:);
        B(2,2:2:8) = dNxy(2,:);
        B(3,1:2:7) = dNxy(2,:);
        B(3,2:2:8) = dNxy(1,:);
        Ke = Ke + B'*C*B*t*det(J);
    end
    dofs = reshape(dof(elem(e,:),:)',[],1);
    Kglobal(dofs,dofs) = Kglobal(dofs,dofs) + Ke;
end
Kglobal = sparse(Kglobal);

%% BC
fijo = [1 1
        0 0
        1 0
        0 0];
fijo = reshape(logical(fijo'),[],1);
libre = ~fijo;

%% Cargas
P = [0 0
     1 0
     0 0
     1 0]*1000; %[N]
P = reshape(P',[],1);
%% Reducción
Kglobalred = Kglobal(libre,libre);
Pred = P(libre);
%% Solver
Dred = Kglobalred\Pred;
D = zeros(doftot,1);
D(libre) = Dred;

P = Kglobal*D;

























