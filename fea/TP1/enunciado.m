
%Types
%1=barra
%2=viga comun
%3=viga con bisagra al comienzo
%4=viga con bisagra al final
apoyos_simples=[7 8];
nod=[0 0;20 0; 20 -8;0 -16;15 -(15+8);%nodos 1 2 3 4 5
    (20+16) -(8+15+12);(16+20+24) -(6+12+15+8);(20+16+24+6) -(8+15+12+6+24)];
elenod=[1 2;1 4;3 4;2 3;3 5;4 5;5 6;6 7;6 8];
eletype=[1 1 1 3 4 2 2 2 1];
ndof=3;
vigasinteresantes=6:9;
elementos=genelementos(elenod);
[Ne,~]=size(elenod);
[N,~]=size(nod);
Ndof=N*ndof;
CB=false(Ndof,1);
%Fuerzas
P=3200;
Fa=-P*15/16;
Fb=-Fa;
R=zeros(Ndof,1);
R(elementos(1,1:2))=[Fa -P/2];
R(elementos(2,4:5))=[Fb -P/2];
Sye=50e3*ones(Ne,1); %limite de fluencia en psi
%Datos Materiales
Ee=30e6*ones(Ne,1);
Ae=8*ones(Ne,1);
ce=4*ones(Ne,1);
he=ce*2;
be=1*ones(Ne,1);
Ie=800*ones(Ne,1);
%
Ie(2)=Ie(2)*4;
Ie(6)=4.890460502746170e-01;
Ie(7)=9.953628668933730e+00;
Ie(8)=1.526000776376096e-01;