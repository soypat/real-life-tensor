%% p2014T3
%% EJ 1
%a. Dependen de los factores geométricos del como el area, la matriz de
%inercia y el largo. También de las propiedades del material, como el
%modulo de young y el modulo de corte transversal.
%b. El vector de cargas depende de las fuerzas tanto internas como externas
%que modifican los desplazamientos.
%% EJ 2
clear; close all; clc
%% Datos
b = 100;
h = 50;
A = b*h;
I = b*h^3/12;
E = 2E5;
M = 2E6;
q = 2;
%% Nodos y elementos
cnod = [0 0
        1 0
        4 0
        7 4
        8 4]*1000;
nnod = 5;
elem = [1 2
        2 3
        3 4
        4 5];
nelem = 4;
%% DOF
dofpornodo = 3;
dof = reshape((1:1:dofpornodo*nnod)',dofpornodo,nnod)';
%% Matriz global
Kglobal = zeros(dofpornodo*nnod);
Klocal = zeros(6);
memlong = zeros(nelem,1);
memT = zeros(6,6,nelem);
for e = 1:nelem
    v = cnod(elem(e,2),:) - cnod(elem(e,1),:);
    long = norm(v);
    memlong(e) = long;
    vd = v/long;
    lambda = [vd 0; -vd(2) vd(1) 0; 0 0 1];
    T = blkdiag(lambda,lambda);
    memT(:,:,e) = T;
    
    X=A*E/long;
    Y4=2*E*I/long;
    Y3=Y4*2;
    Y2=Y4*3/long;
    Y1=Y2*2/long;
    Kbarra=X*[1 -1; -1 1];
    Kviga=[Y1 Y2 -Y1 Y2
           Y2 Y3 -Y2 Y4
          -Y1 -Y2 Y1 -Y2
           Y2 Y4 -Y2 Y3];
    Klocal([1 4],[1 4]) = Kbarra;
    Klocal([2 3 5 6],[2 3 5 6]) = Kviga;
    Kelemento = T'*Klocal*T;
    dofr = [dof(elem(e,1),:) dof(elem(e,2),:)];
    Kglobal(dofr,dofr) = Kglobal(dofr,dofr)+Kelemento;
end
%% BC
fijomatriz = [1 1 1
              1 1 0
              0 0 0
              0 0 0
              1 1 1];
fijo = 1&reshape(fijomatriz',[],1);
libre = ~fijo;
%% Cargas
Pmatriz = [0 0 0
           0 0 -M
           0 -q*memlong(3)/2 -q*memlong(3)^2/12
           0 -q*memlong(3)/2 q*memlong(3)^2/12
           0 0 0];
P = reshape(Pmatriz',[],1);
%% Reducción
Kglobalred = Kglobal(libre,libre);
Pred = P(libre);
%% Solver
Dred = Kglobalred\Pred;
D = zeros(dofpornodo*nnod,1);
D(libre) = Dred;
Dmatriz = reshape(D,dofpornodo,nnod)'
%% Desplazamientos
Delem3 = [Dmatriz(elem(3,1),:) Dmatriz(elem(3,2),:)];
Delem3local = Delem3*memT(:,:,3);
n = 10;
X = linspace(0,memlong(3),n);
Desp3 = zeros(3,n);
% Axiales
N1a =@(x) x/memlong(3);
Na = [-N1a(X); N1a(X)];
Desp3(1,:) = Delem3local([1 4])*Na;
% Transversales
N1t=@(x) 1-3*(x/memlong(3)).^2+2*(x/memlong(3)).^3;
N2t=@(x) x-2*x.^2/memlong(3)+x.^3/memlong(3)^2;
N3t=@(x) 3*(x/memlong(3)).^2-2*(x/memlong(3)).^3;
N4t=@(x) -x.^2/memlong(3) +x.^3/memlong(3)^2;
Nt = [N1t(X); N2t(X); N3t(X); N4t(X)];
Desp3(2,:) = Delem3local([2 3 5 6])*Nt;
% Angular
Desp3(3,:) = ones(1,n);
Desp3rotado = (memT(1:3,1:3,3)'*Desp3)'
%% Tensiones en A
Delem4 = [Dmatriz(elem(4,1),:) Dmatriz(elem(4,2),:)];
X = 0;
N1=@(x) -6/memlong(4)^2+12*x/memlong(4)^3;
N2=@(x) -4/memlong(4)+6*x/memlong(4)^2;
N3=@(x) 6/memlong(4)^2-12*x/memlong(4)^3;
N4=@(x) -2/memlong(4)+6*x/memlong(4)^2;
B = [N1(X); N2(X); N3(X); N4(X)];
SigmaFlex = h/2*E*Delem4([2 3 5 6])*B
