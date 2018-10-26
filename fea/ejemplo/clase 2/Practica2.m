%% Practica 2
% Datos geometricos
clear all; close all; clc
L=20; A=10000; E=210000; F=2500;
c=cos(pi()/6); s=sin(pi()/6);
%% Nodos
nod=[1 2 3 4];
cnod=[0 0 0
    -20 0 0
    -10 20 20*c
    -10 0 20*c];
elem=[1 2
    1 3
    1 4
    2 3
    2 4
    3 4];
%% Matriz global
Kg=zeros(12,12);
dof=[1 2 3; 4 5 6; 7 8 9; 10 11 12];
for i=1:6
    v=cnod(elem(i,2),:)-cnod(elem(i,1),:);
    Ke=rotador(v,A,E);
    dofs = [dof(elem(i,1),:),dof(elem(i,2),:)];
    Kg(dofs,dofs) = Kg(dofs,dofs) + Ke;
end
R=[0 0 -2500].';
free=dof(3,:);

free=1&[0 0 0 0 0 0 1 1 1 0 0 0];
Kgr=Kg(free,free);
D=Kgr\R;

