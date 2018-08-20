
%Types
%1=barra
%2=viga comun
%3=viga con bisagra al comienzo
%4=viga con bisagra al final
apoyos_simples=[7 8];
nod=[0 0;20 0; 20 -8;0 -16;15 -(15+8);%nodos 1 2 3 4 5
    (20+16) -(8+15+12);(16+20+24) -(6+12+15+8);(20+16+24+6) -(8+15+12+6+24)];
elenod=[1 2;1 4;3 4;2 3;3 5;4 5;5 6;6 7;6 8];
eletype=[1 1 1 3 4 3 2 2 1];
ndof=3;
vigasinteresantes=[6:8];
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
Sye=77e3*ones(Ne,1); %limite de fluencia en psi

%Datos Materiales
safetyfactor=8;
Ie=zeros(Ne,1);
Ee=30e6*ones(Ne,1);
Ee(2)=Ee(2)*6;
Ae=zeros(Ne,1);

% he=[1 4 1 4 4 ...
%     4 4 6 4];
% 
% be=[1 2 .5 1 2 ...
%     1.5 2 2 1];%Vigas 6 7 8 9


% he=[1 4 1 4 4 ...
%     4 5 1 4];
% 
% be=[.75 2 .5 1 2 ...
%     1.5 1.5 1.5 1.5];%Vigas 6 7 8 9
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
