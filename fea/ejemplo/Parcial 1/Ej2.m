%% Ej2
clear; close all; clc
%% Datos
a = 40;
h = 120;
A = a*h; 
I = a*h^3/12;
E = 210000;
alfa = 1E-6;
tita = 100;
F = -750;
%% Nodos y elementos
cnod = [0 0
        0 280
        112.58 345
        112.58+25.98 360
        199.19 396];
nnod = 5;
elem = [1 2
        2 3
        3 4
        4 5];
nelem = 4;
%% Dof
dofpornodo = 3;
dof = reshape((1:1:dofpornodo*nnod)',dofpornodo,nnod)';
%% Matriz de rigidez
Kparcial = zeros(dofpornodo*nnod);
Kviga2D = zeros(6);
memT = zeros(6,6,nelem);
memlong = zeros(nelem,1);
for e = 2:nelem
    v = cnod(elem(e,2),:) - cnod(elem(e,1),:);
    long = norm(v);
    memlong(e) = long;
    vd = v/long;
    lambda = [vd 0; -vd(2) vd(1) 0; 0 0 1];
    T = blkdiag(lambda,lambda);
    memT(:,:,e) = T;
    X = E*A/long;
    Y4 = 2*E*I/long;
    Y3 = Y4*2;
    Y2 = Y4*3/long;
    Y1 = Y2*2/long;
    Kbarra = X*[1 -1;-1 1];
    Kviga = [Y1 Y2 -Y1 Y2
             Y2 Y3 -Y2 Y4
            -Y1 -Y2 Y1 -Y2
             Y2 Y4 -Y2 Y3];
    Klocal([1 4],[1 4]) = Kbarra;
    Klocal([2 3 5 6],[2 3 5 6]) = Kviga;
    Kelem = T'*Klocal*T;
    dofr = [dof(elem(e,1),:) dof(elem(e,2),:)];
    Kparcial(dofr,dofr) = Kparcial(dofr,dofr)+Kelem;
end
Kparcial([2,5],[2,5]) = Kparcial([2,5],[2,5]) + E*A/280*[1 -1;-1 1];
Kglobal = Kparcial([2 4:15],[2 4:15]);
%% BC
fijo = [1 zeros(1,9) 1 1 1]';
libre = ~fijo;
%% Carga
P = zeros(13,1);
P([6 9]) = F*[1; 1];
FT = alfa*tita*E*A;
P(3) = FT;
%% Reducción
Kglobalred = Kglobal(libre,libre);
Pred = P(libre);
%% Solver
Dred = Kglobalred\Pred;
D = [0 Dred' 0 0 0]';
%% Tensiones
Dpornodo = reshape(D(2:end),dofpornodo,nnod-1)';
SigmaAxial = zeros(nelem-1,1);
SigmaFlexSup = zeros(nelem-1,30);
SigmaTotalSup = SigmaFlexSup;
SigmaTotalInf = SigmaFlexSup;
memsub = SigmaFlexSup;
for e = 1:nelem-1
    Delemento = [Dpornodo(elem(e,1),:) Dpornodo(elem(e,2),:)]';
    Dlocal = memT(:,:,e+1)*Delemento;
    Ba = [-1 1]/memlong(e+1);
    SigmaAxial(e) = E*Ba*Dlocal([1 4]);   
    sub = 0:memlong(e+1)/29:memlong(e+1);
    memsub(e,:) = sub;
    N1=@(x) -6/memlong(e+1)^2+12*x/memlong(e+1)^3;
    N2=@(x) -4/memlong(e+1)+6*x/memlong(e+1)^2;
    N3=@(x) 6/memlong(e+1)^2-12*x/memlong(e+1)^3;
    N4=@(x) -2/memlong(e+1)+6*x/memlong(e+1)^2;
    Bf = [N1(sub); N2(sub); N3(sub); N4(sub)]';
    SigmaFlexSup(e,:) = (h/2*E*Bf*Dlocal([2 3 5 6]))';
    SigmaFlexInf = -SigmaFlexSup;
    SigmaTotalSup(e,:) = SigmaFlexSup(e,:)+SigmaAxial(e);
    SigmaTotalInf(e,:) = SigmaFlexInf(e,:)+SigmaAxial(e);
end