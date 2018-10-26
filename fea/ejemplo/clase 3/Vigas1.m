%% Vigas sin axial 2D
clear all; close all; clc
%% Datos
E=207000;
a=20;
Iz=(a^4)/12;
L=1000;
F=300;
%% Nodos
cnod=[0 0
    L 0];
nnod=size(cnod,1);
elem=[1 2];
nelem=size(elem,1);
%% Dof
dofpornodo=2; %Grados de libertad por nodo
dof=(1:1:dofpornodo*nnod)';
dof=(reshape(dof,dofpornodo,nnod))';
%% Matriz global
Kg=zeros(dofpornodo*nnod);
for e=1:nelem
    long=norm(cnod(elem(e,2),:)-cnod(elem(e,1),:));
    k1=2*E*Iz/long;
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
libre=~fijo;
libre=reshape(libre.',[],1);
%% Cargas
P=[0 0 -F 0];
%% Reducción de matriz global
Kgr=Kg(libre,libre);
Pr=P(libre)';
Dr=Kgr\Pr;
