%example58  in, lbf,  psi, etc
Iz=100;%in^4
Iy=100;
K=50;
E=30e3;
nu=.25;
L=100;
A=10;
nodos=[];
index=[1 2 3 4 5 6 7 8 9 10 11 12;
    13 14 15 16 17 18, 7 8 9 10 11 12;
    19 20 21 22 23 24, 7 8 9 10 11 12];
klocal=Kvuw(E,A,Iz,Iy,K,nu,L);

T1=Tvuw([1 0 0],[0 1 0]);
T2=Tvuw([0 1 0],[1 0 0]);
T3=Tvuw([0 0 1],[0 1 0]);
krot1=T1'*klocal*T1;
krot2=T2'*klocal*T2;
krot3=T3'*klocal*T3;

Ndof=6*4;
kG=zeros(Ndof);
kG(index(1,:),index(1,:))=kG(index(1,:),index(1,:))+krot1;
i=2;
kG(index(i,:),index(i,:))=kG(index(i,:),index(i,:))+krot2;
i=3;
kG(index(i,:),index(i,:))=kG(index(i,:),index(i,:))+krot3;

% Ktotal=krot1+krot2+krot3;
CB=false(24,1);
nodo=@(n) [n*6-5 n*6-4 n*6-3 n*6-2 n*6-1 n*6];
CB(nodo(1))=true;
CB(nodo(3))=true;
CB(nodo(4))=true;
R=zeros(6,1);
% R(2)=-50;
XX=false(Ndof,1);
XX(8)=true;
CC=~XX;
XXr=XX(~CB);
R(4)=1e3;
Kr=kG(~CB,~CB);
U=Kr\R;
D=zeros(24,1);
D(~CB)=U;

