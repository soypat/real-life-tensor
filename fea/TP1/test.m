clear all
%Types
%1=barra
%2=viga comun
%3=viga con bisagra al comienzo
%4=viga con bisagra al final
ndof=3;

nod=[0 0;1 1;2 2];
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
Sye=50e3*ones(Ne,1); %limite de fluencia en psi
vigasinteresantes=[1,2];

safetyfactor=2;
he=8*ones(Ne,1);
be=1*ones(Ne,1);
% for i=1:Ne
%     if eletype(i)==1
%         Ae(i)=be(i)^2*pi/4;
%         Ie(i)=be(i)^4*pi/(64);
%     else
%         Ae(i)=be(i)*he(i);
%         Ie(i)=(be(i)*he(i)^3)/12;
%     end
% end
Ie=11.3*ones(Ne,1);
Ae=3.83*ones(Ne,1);

