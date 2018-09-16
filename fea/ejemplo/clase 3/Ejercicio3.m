%% Ejercicio ultimo
clear all; close all; clc
%% Datos
E=210000;
I1=700E-6*1000^4;
I=[I1 2*I1];
L1=12000;
L2=6000;
F=-100000;
q=-10;
%% Nodos
nelem1=10;
nelem2=5;
cnod=[0:L1/nelem1:L1 L1+L2/nelem2:L2/nelem2:L1+L2];
nnod=size(cnod,2);
elem=[1:1:nelem1+nelem2; 2:1:nelem1+nelem2+1]';
%% Dof
dofpornodo=2; %Grados de libertad por nodo
dof=(1:1:dofpornodo*nnod)';
dof=(reshape(dof,dofpornodo,nnod))';
%% Matriz global
Kg=zeros(dofpornodo*nnod);
for e=1:nelem1
    long=cnod(elem(e,2))-cnod(elem(e,1));
    k1=2*E*I(1)/long;
    k2=k1*2;
    k3=k1*3/long;
    k4=k3*2/long;
    Ke=[k4 k3 -k4 k3
        k3 k2 -k3 k1
        -k4 -k3 k4 -k3
        k3 k1 -k3 k2];
    dofr=[dof(elem(e,1),:) dof(elem(e,2),:)];
    Kg(dofr,dofr)=Kg(dofr,dofr)+Ke;
end
for e=nelem1+1:nelem1+nelem2
    long=cnod(elem(e,2))-cnod(elem(e,1));
    k1=2*E*I(2)/long;
    k2=k1*2;
    k3=k1*3/long;
    k4=k3*2/long;
    Ke=[k4 k3 -k4 k3
        k3 k2 -k3 k1
        -k4 -k3 k4 -k3
        k3 k1 -k3 k2];
    dofr=[dof(elem(e,1),:) dof(elem(e,2),:)];
    Kg(dofr,dofr)=Kg(dofr,dofr)+Ke;
end
%% BC
fijo=false(nnod,dofpornodo);
fijo(1,1:2)=true;
fijo(nelem1+1,1)=true;
libre=~fijo;
libre=reshape(libre.',[],1);
%% Cargas
P=zeros(nnod,dofpornodo);
for e=1:nelem1
    long=cnod(elem(e,2))-cnod(elem(e,1));
    P(e,1)=P(e,1)+q*long/2;
    P(e+1,1)=P(e+1,1)+q*long/2;
end
P(nelem1+1,2)=-q*(long^2)/12;
P(end,1)=F;
P=reshape(P',[],1);
%% Reducción de matriz global
Kgr=Kg(libre,libre);
Pr=P(libre);
%% Solver
Dr=Kgr\Pr;
D=zeros(nnod*2,1);
D(libre)=Dr;
d=reshape(D,2,nnod);
plot(cnod,d(1,:))
grid on
%% Reacciones
R=Kg*D;
%% Funciones de forma
N1=@(x) 1-3*(x/long).^2+2*(x/long).^3;
N2=@(x) x-2*x.^2/long+x.^3/long^2;
N3=@(x) 3*(x/long).^2-2*(x/long).^3;
N4=@(x) -x.^2/long +x.^3/long^2;

