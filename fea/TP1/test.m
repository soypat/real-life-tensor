clear all
%Types
%1=barra
%2=viga comun
%3=viga con bisagra al comienzo
%4=viga con bisagra al final
ndof=3;

nod=[0 0;1 1;0 2];
elenod=[1 2;2 3];
eletype=[2 2];
apoyos_simples=[1 3];
elementos=genelementos(elenod);

[Ne,~]=size(elenod);
[N,~]=size(nod);
Ndof=N*ndof;

%Fuerzas
P=3200;
% Fa=-P*15/16;
% Fb=-Fa;
R=zeros(Ndof,1);
R=fuerzapuntual(R,2,1e4,-1e4,0);
% R(elementos(1,1:2))=[Fa -P/2];
% R(elementos(2,4:5))=[Fb -P/2];

%Datos Materiales
Ee=30e6*ones(Ne,1);%Modulo young psi
Ae=8*ones(Ne,1);
ce=4*ones(Ne,1);
be=4*ones(Ne,1);
Ie=800*ones(Ne,1);
Sye=50e3*ones(Ne,1); %limite de fluencia en psi
vigasinteresantes=[1,2];
