L=1;
alto=sqrt(2)/(2*L);
nod=[0 0;alto alto;L+alto alto;alto+2*L alto;L 0;L+alto/2 alto/2;2*L 0;L+alto alto];
elenod=[1 2;2 3;3 4;5 6;6 8;7 6];
[Ne,~]=size(elenod);
Le=zeros(Ne,1);
for i = 1:Ne
    nodestart=elenod(i,1);
    nodeend=elenod(i,2);
    lx=nod(nodeend,1)-nod(nodestart,1);
    ly=nod(nodeend,2)-nod(nodestart,2);
    Le(i)=sqrt(lx^2+ly^2);
    phide(i)=atan2d(ly,lx);%angulo en degrees
end
nodeDofs=[1 2 0;
    3 4 5;
    6 7 8;
    9 10 11;
    12 13 14;
    15 16 17;
    18 19 0;
    6 7 20];
Ndof=max(max(nodeDofs));
E=200e9;
b=0.04;%m
A=b^2;
Ap=b^2*pi/4;
I=b^4/12;
kG=zeros(Ndof);

i=1;
kbarra=Kb(E,A,Le(i));
T=Tb(phide(i));
kbarrarotada=T*kbarra*T';%LAS BARRAS QUEDAN tamaño 4x4
elementos=[nodeDofs(elenod(i,1),[1 2]) nodeDofs(elenod(i,2),[1 2])];
kG(elementos,elementos)=kG(elementos,elementos)+kbarrarotada;

i=6;
barnod=[1 2 4 5];
kbarra=Kb(E,Ap,Le(i));
T=Tb(phide(i));
kbarrarotada=T*kbarra*T';
elementos=[nodeDofs(elenod(i,1),[1 2]) nodeDofs(elenod(i,2),[1 2])];
kG(elementos,elementos)=kG(elementos,elementos)+kbarrarotada;
for i=[2 3 4 5]
    klocal=Kv(E,A,I,Le(i));
    T=Tvu(phide(i));
    klocalrotada=T'*klocal*T;
    elementos=[nodeDofs(elenod(i,1),:) nodeDofs(elenod(i,2),:)];
    kG(elementos,elementos)=kG(elementos,elementos)+klocalrotada;
end
CB=false(Ndof,1);
CB([1 2 12 13 18 19])=true;
Kr=kG(~CB,~CB);
R=zeros(Ndof,1);
R(10)=-1000;%N
F=R(~CB);
U=Kr\F
function [T] = Tb(nine)
T=[cosd(nine) 0;
    sind(nine) 0;
    0 cosd(nine);
    0 sind(nine)];
end
