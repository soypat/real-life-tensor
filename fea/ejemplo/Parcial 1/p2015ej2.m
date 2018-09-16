%% Parcial 2015 ej 2
clear; close all; clc
%% Datos
b = 50; %mm
h = 100; %mm
A = b*h; %mm^2
I = b*h^3/12; %mm^4
E = 70000; %MPa
alfa = 1E-6; 
tita = 100;
L = 2000; %mm
%% Nodos
cnod = [0 0
        0 -L
        L -L];
nnod = 3;
elem = [1 2
        2 3];
nelem = 2;
%% Dof
dofpornodo = 3;
dof = [1 2 3
       4 5 6
       7 8 9];
%% Matriz global
Kglobal = zeros(dofpornodo*nnod);
Klocal = zeros(6);
X = A*E/L;
Y4 = 2*E*I/L;
Y3 = Y4*2;
Y2 = Y4*3/L;
Y1 = Y2*2/L;
Kbarra = X*[1 -1;-1 1];
Kviga = [Y1 Y2 -Y1 Y2
         Y2 Y3 -Y2 Y4
        -Y1 -Y2 Y1 -Y2
         Y2 Y4 -Y2 Y3];
Klocal([1 4],[1 4]) = Kbarra;
Klocal([2 3 5 6],[2 3 5 6]) = Kviga;
lambda1 = [0 -1 0; 1 0 0; 0 0 1];
T1 = blkdiag(lambda1,lambda1);
K1 = T1'*Klocal*T1;
K2 = Klocal;
Kglobal(1:6,1:6) = Kglobal(1:6,1:6)+K1;
Kglobal(4:9,4:9) = Kglobal(4:9,4:9)+K2;
%% BC
fijo = [1 1 1 0 0 0 1 1 1]';
libre = ~fijo;
%% Cargas
FT = alfa*tita*E*A;
P = [0 0 0 -FT 0 0 FT 0 0]';
%% Reducción
Kglobalred = Kglobal(libre,libre);
Pred = P(libre);
%% Solver
Dred = Kglobalred\Pred;
D = [0 0 0 Dred' 0 0 0]';
%% Reacciones
Pp = Kglobal*D-P;
%% Tensiones
D1 = D(1:6);
D1 = T1*D1;
D2 = D(4:9);
% Axiales
Baxial = [-1 1]/L;
SigmaAxial1 = E*Baxial*D1([1 4]);
SigmaAxial2 = E*Baxial*D2([1 4])-FT/A;
% Flexión (en el trocen)
N1 = 0;
N2 = -1/L;
N3 = 0;
N4 = 1/L;
Bflex = [N1 N2 N3 N4];
SigmaFlex1Sup = -h/2*E*Bflex*D1([2 3 5 6]);
SigmaFlex2Sup = -h/2*E*Bflex*D2([2 3 5 6]);
SigmaFlex1Inf = -SigmaFlex1Sup;
SigmaFlex2Inf = -SigmaFlex2Sup;
Sigma1Sup = SigmaFlex1Sup+SigmaAxial1;
Sigma2Sup = SigmaFlex2Sup+SigmaAxial2;
Sigma1Inf = SigmaFlex1Inf+SigmaAxial1;
Sigma2Inf = SigmaFlex2Inf+SigmaAxial2;
SigmaViga1 = [Sigma1Sup SigmaAxial1 Sigma1Inf]'
SigmaViga2 = [Sigma2Sup SigmaAxial2 Sigma2Inf]'
%% Deformación
X = linspace(0,L,30);
N1=@(x) 1-3*(x/L).^2+2*(x/L).^3;
N2=@(x) x-2*x.^2/L+x.^3/L^2;
N3=@(x) 3*(x/L).^2-2*(x/L).^3;
N4=@(x) -x.^2/L +x.^3/L^2;
N = [N1(X); N2(X); N3(X); N4(X)];
Desp1 = D1([2 3 5 6])'*N;
Desp2 = D2([2 3 5 6])'*N;