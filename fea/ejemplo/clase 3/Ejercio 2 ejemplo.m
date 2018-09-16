%% Ejercicio 2, ejemplo
clear all; close all; clc
%% Datos
E=210000;
I=2E-4*1000^4;
L=3000;
F=-50000;
%% Nodos
cnod=[0 0
      0 L
      0 2*L
      -1 2*L];
nnod=size(cnod,2);
elem=[1 2
      2 3
      3 4];
%% Dof
dofpornodo=2; %Grados de libertad por nodo
dof=(1:1:dofpornodo*nnod)';
dof=(reshape(dof,dofpornodo,nnod))';
%% Matriz global
Kg=zeros(dofpornodo*nnod);