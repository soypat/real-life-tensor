%% Vigas con un pin
clear; close all; clc
%% Datos
b = 100;
h = 100;
A = b*h;
I = b*h^3/12;
E = 2E5;
P = 10000;
%% Nodos
cnod = [0 -1
        0 0
        1 0
        2 0]*1000;
nnod = 4;
elem = [1 2
        1 3
        2 3
        3 4];
nelem = 4;
%% DOF
dofpornodo = 3;
dof = reshape((1:1:dofpornodo*nnod)',dofpornodo,nnod)';
%% Matriz de rigidez
Kglobal = zeros(dofpornodo*nnod);
for e = 1:nelem
    v=cnod(elem(e,2),:)-cnod(elem(e,1),:);
    long=norm(v);
    vd=v/long;
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
    lambda=[vd(1) vd(2) 0; -vd(2) vd(1) 0; 0 0 1];
    T=blkdiag(lambda,lambda);
    Kelem=T'*Klocal*T;
    dofr=[dof(elem(e,1),:) dof(elem(e,2),:)];
    Kglobal(dofr,dofr)=Kglobal(dofr,dofr)+Kelem;
end
%% Pin
% Existe un pin en el nodo 3 que une el elemento 2 con el 3 y el 4
Kglobal(dof(elem(2,2),3),dof(elem(2,2),3)) = 0;