L=1000;%mm
E=210e3;%MPa
a=20;%mm

% alto=sqrt(2)/(2*L);
nod=[0 0;0 L;0 2*L;-L L;-L 2*L;-1.5*L 0;-1.5*L L;-1.5*L 2*L];
N=length(nod);
R=zeros(length(nod)*3,1); %vector fuerzas, 
elenod=[8 5;5 3;7 4;4 2;7 6;5 4;5 2;4 1;3 2;2 1];
CB=false(N*3,1);
apoyos_simples=[];
empotramientos=[1 6];
eletype=ones(1,length(elenod))*2;
eletype(5)=22;%simetría, 
%CASO 1
P=10e3;
% R=fuerzapuntual(R,4,0,-1000,0);%Aplico 1000N sobre mi nodo 4 en direccion y, sentido negativo
%es lo mismo que R(11)=-1000;%N
nodo=@(n) [n*3-2 n*3-1 n*3];
%CASO 1 SIMETRICO
R=fuerzapuntual(R,3,P/2,0,0);
R=fuerzapuntual(R,8,0,P/2,0);
CB(nodo(8))=~~[1 0 1]';
CB(nodo(7))=~~[1 0 1]';

% CASO 2 ANTISIMETRICO, LO EVALUO AL FINAL
% R=fuerzapuntual(R,3,P/2,0,0);
% CB(nodo(8))=~~[0 1 0]';
% CB(nodo(7))=~~[0 1 0]';





A=a^2;
I=a^4/12;


[Ne,~]=size(elenod);%Obtengo cantidad de elementos usados Ne



safetyfactor=2;
Ee=E*ones(Ne,1);
Ae=A*ones(Ne,1);
Ie=I*ones(Ne,1);

he=a*ones(Ne,1);
be=a*ones(Ne,1);
Sye=350e6*ones(Ne,1);

% Ie(5)=I/2;
% Ae(5)=A/2;

graficar=true;
showresults=false;
% graficar=false;
vigasinteresantes=[];

G2DFEAS
Dc1=D;%Guardo desplazamientos de simetría.
R=fuerzapuntual(R,3,P/2,0,0);
CB(nodo(8))=~~[0 1 0]';
CB(nodo(7))=~~[0 1 0]';

G2DFEAS
Dc2=D;%Guardo desplazamientos de antisimetría

Dt=Dc1+Dc2;
disp(Dt(nodo(3)))