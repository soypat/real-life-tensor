%% EJ1
clear; close all; clc
tic
%% Datos
% Tensión plana
E = 30E6; % [psi]
nu = 0.25;
t = 0.5; % [in]
p = 300; % [psi]
f = 1000; % [lb]
% 7 8 9
% 4 5 6
% 1 2 3
%% Nodos
cnod = [0	0
        1.5	0
        3	0
        0	1
        1.5	1
        3	1
        0	2
        1.5	2
        3	2]; % [in]
nnod = size(cnod,1);
elem = [1	2	5	4
        2	3	6	5
        4	5	8	7
        5	6	9	8];
elempornodo = size(elem,2);
nelem = size(elem,1);
%% DOF
dofpornodo = 2;
doftot = dofpornodo*nnod;
dof = reshape((1:1:doftot)',dofpornodo,nnod)';
%% BC
fijos = false(nnod,dofpornodo);
% fijos([1 4 7],:) = true;
fijos([1 7],:) = true;
fijos(3,2) = true;
% fijos(2,2) = true;
fijos = reshape(fijos',[],1);
libre = ~fijos;
%% Cargas
P = zeros(nnod,dofpornodo);
P(7,2) = -p*.75*t;
P(8,2) = -p*.75*t*2;
P(9,2) = -p*.75*t-f;
P = reshape(P',[],1);
%% Matriz
% Matriz constitutiva
C = [1 nu 0
     nu 1 0
     0 0 (1-nu)/2]*E/(1-nu^2)*t;
% Funciones de forma
syms as bs x y
N1 = (as-x)*(bs-y)/4/as/bs;
N2 = (as+x)*(bs-y)/4/as/bs;
N3 = (as+x)*(bs+y)/4/as/bs;
N4 = (as-x)*(bs+y)/4/as/bs;
N = [N1 N2 N3 N4];
% Derivadas respecto a x e y
dNx = diff(N,x);
dNy = diff(N,y);
% Matriz de deformación
B = [dNx(1)   0    dNx(2)   0    dNx(3)   0    dNx(4)   0
       0    dNy(1)   0    dNy(2)   0    dNy(3)   0    dNy(4)
     dNy(1) dNx(1) dNy(2) dNx(2) dNy(3) dNx(3) dNy(4) dNx(4)];
Kglobal = zeros(doftot);
for e = 1:nelem
    a = norm(cnod(elem(e,1),:)-cnod(elem(e,2),:))/2;
    b = norm(cnod(elem(e,2),:)-cnod(elem(e,3),:))/2;
    Be = subs(B,{as,bs},{a,b});
    Ke = int(int(Be'*C*Be,x,-a,a),y,-b,b);
    dofselec = reshape((dof(elem(e,:),:))',[],1);
    Kglobal(dofselec,dofselec) = Kglobal(dofselec,dofselec)+Ke;
end
Kglobalvisual = Kglobal;
Kglobal = sparse(Kglobal);
%% Solver
Dr = Kglobal(libre,libre)\P(libre);
D = zeros(doftot,1);
D(libre) = Dr;
Dvisual = reshape(D,dofpornodo,nnod)'
%% Graficos
cnodfinal = cnod+Dvisual*1000;
figure
hold on; grid on; axis equal;
plot(cnod(:,1),cnod(:,2),'.')
plot(cnodfinal(:,1),cnodfinal(:,2),'.')
toc