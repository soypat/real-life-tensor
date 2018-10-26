%% Ejercicio Pórtico 2D
clear all ; close all ; clc
%% Datos
d2x=5.7; %mm
d2y=-0.0224; %mm
phi2=0.00523; %rad
L1= 4000; %mm
L2= 10000; %mm
q= 300; %N/mm
E= 210000; %MPa
I= 0.5E-4*1000^4; %mm^4
A= 0.5E-2*1000^2; %mm^2
%% Nodos
cnod=[0 0
      0 L1
      L2 L1
      0 L1*2
      L2 L1*2];
nnod=size(cnod,1);
elem=[1 2
      2 3
      2 4
      4 5];
nelem=size(elem,1);
%% DOF
dofpornodo=3;
dof=reshape([1:1:dofpornodo*nnod]',dofpornodo,nnod)';
%% Matriz global
Kg=zeros(nnod*dofpornodo,nnod*dofpornodo);
Kl=zeros(6,6);
for e=1:nelem
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
    Kl([1 4],[1 4])=Kbarra;
    Kl([2 3 5 6],[2 3 5 6])=Kviga;
    dofr=[dof(elem(e,1),:) dof(elem(e,2),:)];
    lambda=[vd(1) vd(2) 0; -vd(2) vd(1) 0; 0 0 1];
    T=blkdiag(lambda,lambda);
    Ke=T'*Kl*T;
    Kg(dofr,dofr)=Kg(dofr,dofr)+Ke;
end
%% BC
fijo=zeros(nnod,dofpornodo);
fijo([1 3 5],:)=1;
libre=~fijo;
libre=reshape(libre.',[],1);
%% Cargas
P=zeros(nnod,dofpornodo);
izqtriangulo=[7*L1/20 0 -L1^2/20];
dertriangulo=[3*L1/20 0 L1^2/30];
izqcuadrado=[L1/2 0 -L1^2/12];
dercuadrado=[L1/2 0 L1^2/12];
P([1 2 4],:)=q/2*[izqtriangulo+izqcuadrado
              dertriangulo+dercuadrado+izqtriangulo
              dertriangulo];
P=reshape(P',[],1);
%% Reductor
Kgr=Kg(libre,libre);
Pr=P(libre);
Dr=Kgr\Pr
%% Desplazamientos desconocidos
% % Los conocidos son todos los dof de el nodo 2
% doffila=reshape(dof',[],1);
% despcono=1&[ones(1,9) 0 0 0 ones(1,3)];
% despdesc=~despcono;
% Kgxc=Kg(doffila(despcono),doffila(despdesc));
% P
