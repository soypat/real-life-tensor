L=1;
alto=sqrt(2)/(2*L);
nod=[0 0;alto alto;L+alto alto;alto+2*L alto;L 0;alto/2+L alto/2;L*2 0];
R=zeros(length(nod)*3,1); %vector fuerzas, 
elenod=[1 2;2 3;3 4;5 6;6 3;7 6];

eletype=[1 2 2 2 4 1];
apoyos_simples=[1 5 7];

R=fuerzapuntual(R,4,0,-1000,0);%Aplico 1000N sobre mi nodo 4 en direccion y, sentido negativo
%es lo mismo que R(11)=-1000;%N

empotramientos=[];

CB=false(22,1);

E=200e9;%Pa
a=.04;%m
A=a^2;
I=a^4/12;
Ip=a^4*pi/64;
Ap=a^2*pi/4;

[Ne,~]=size(elenod);%Obtengo cantidad de elementos usados Ne

safetyfactor=2;
Ee=E*ones(Ne,1);
Ae=A*ones(Ne,1);
Ie=I*ones(Ne,1);
Ae(6)=Ap;
he=a*ones(Ne,1);
be=a*ones(Ne,1);
Sye=350e6*ones(Ne,1);

graficar=true;
vigasinteresantes=[2,3];

G2DFEAS