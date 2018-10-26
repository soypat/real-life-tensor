%Parcial
nod=[0 0;0 .2;.2 .2+.15;.4 .2+.15+.15;.2+.2+.1 .2+.15+.15+.15/2;.2+.2+.1, .2+.15+.15/2;
    .2+.15 .2;.4 -.1;.2 .2];
N=size(nod,1);
fnod=@(n) [n*3-2 n*3-1 n*3];
fnodb=@(n) [n*3-2 n*3-1];
CB=false(27,1);
CB(1)=true;
% scatter(nod(:,1),nod(:,2))
enod=[1 2;2 3;3 4;4 5;4 6;6 5;7 3;7 9;9 3;8 9];
elenod=enod;
Ne=size(enod,1);
eletype=[1 2 2 2 2 1, 4 2 1, 1];% 1 barra, 2 viga, 3 rotula final, 4 piston
apoyos_simples=[7];
empotramientos=[1 8];
Ndof=N*3+1;
R=zeros(Ndof);
P=(-50-100)/2;
M=(8.33333+12.5)/2;
R=fuerzapuntual(R,5,0,P,M);
R=fuerzapuntual(R,6,0,P,M);
Ee=200e9*ones(Ne,1);
b=.05;
h=.15;
be=b*ones(Ne,1);
he=h*ones(Ne,1);
Ae=b*h*ones(Ne,1);
Ae(end)=.05^2*pi/4;
Sye=ones(Ne,1);

Ie=b*h^3/12*ones(Ne,1);
safetyfactor=1;
show_results=false;

G2DFEAS