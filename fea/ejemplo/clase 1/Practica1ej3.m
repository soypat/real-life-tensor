%% Problema 3
clear all; close all; clc
%% Datos
E=210000;
A=500;
L=5000;
F=100000;
k=4;
alfa=pi()/3;
%% Sin simetría
cnod=L*[0 0
    -cos(alfa) sin(alfa)
    cos(alfa) sin(alfa)
    0 -1];
nnod=size(cnod,1);
elem=[1 2
    1 3
    1 4];
nelem=size(elem,1);
dof=[1 2
    3 4
    5 6
    7 8];
Kg=zeros(2*nnod);
%Elementos
for e=1:nelem-1
    v=cnod(elem(e,2),:)-cnod(elem(e,1),:);
    long=norm(v);
    vd=v/long;
    T=[vd 0 0; 0 0 vd];
    Ke=E*A/long*[1 -1; -1 1];
    K=T.'*Ke*T;
    dofr=[dof(elem(e,1),:) dof(elem(e,2),:)];
    Kg(dofr,dofr)=Kg(dofr,dofr)+K;
end
%Resorte
Ke=k*[1 -1; -1 1];
T=[0 -1 0 0; 0 0 0 -1];
K=T.'*Ke*T;
Kg([1 2 7 8],[1 2 7 8])=Kg([1 2 7 8],[1 2 7 8])+K;
Pr=[0 -F].';
Kgr=Kg([1 2],[1 2]);
Dr=Kgr\Pr;
D=[Dr.' 0 0 0 0 0 0].';
P=Kg*D;
disp(strcat(['El desplazamiento vertical sin usar simetria: ',num2str(Dr(2)),' mm']))
%% Con simetría
% Con simetría, se puede suponer que el nodo central esta limitado solo a
% un movimiento vertical
cnod=L*[0 0
    -cos(alfa) sin(alfa)
    0 -1];
nnod=size(cnod,1);
elem=[1 2
    1 3];
nelem=size(elem,1);
dof=[1 2
    3 4
    5 6];
Kg=zeros(2*nnod);
%Elemento
v=cnod(elem(1,2),:)-cnod(elem(1,1),:);
long=norm(v);
vd=v/long;
T=[vd 0 0; 0 0 vd];
Ke=E*A/long*[1 -1; -1 1];
K=T.'*Ke*T;
dofr=[dof(elem(1,1),:) dof(elem(1,2),:)];
Kg(dofr,dofr)=Kg(dofr,dofr)+K;
%Resorte
Ke=k*[1 -1; -1 1];
T=[0 -1 0 0; 0 0 0 -1];
K=T.'*Ke*T;
Kg([1 2 5 6],[1 2 5 6])=Kg([1 2 5 6],[1 2 5 6])+K;
%Solver
Pr=-F/2;
Kgr=Kg(2,2);
Dr=Kgr\Pr;
D=[0 Dr 0 0 0 0].';
P=Kg*D;
disp(strcat(['El desplazamiento vertical usando simetria: ',num2str(Dr),' mm']))