%% Ejercicio 8 el culo te abrocho
format short g
syms x y real
alto=2;
ancho=4;
b=ancho/4; %ancho dividido 4 porque uso 8 elementos
h=alto/2; %Alto dividido 2
a=h;
N1s=(b-x)*(h-y)/(4*b*h); %Eje centrado
N2s=(b+x)*(h-y)/(4*b*h);
N3s= (b+x)*(h+y)/(4*b*h);
N4s= (b-x)*(h+y)/(4*b*h);%Lo saque del Daryl Logan. Está tambien en el Cook
N=[N1s N2s N3s N4s]; %Una fila de matriz de funciones de forma asi puedo derivar split
B=[diff(N,x);diff(N,y)];
C=[25 0;0 25]; %25 W/mk isotropo
k=int(int(B'*C*B,x,-a,a),y,-b,b )
eig(k) % 
%% Comenzamos a buscar la solucion a problema elementos finitos
%Desconocidos
x=[]; %nodos desconocidos
meinDOF=[1 2 6 7;2 3 7 8;3 4 8 9;4 5 9 10;... primera fila de elementos
    6 7 11 12;7 8 12 13;8 9 13 14;9 10 14 15]; %segunda fila de elementos
nodos=[0 0;.25 0;.5 0;.75 0;1 0; ...
       0 .25;.25 .25;.5 .25;.75 .25;1 .25; ...
       0 .5;.25 .5;.5 .5;.75 .5;1 .5];
Ne=8;
Ndof=15;
kG=zeros(Ndof);

for e=1:Ne
    index=meinDOF(e,:);
    kG(index,index)=kG(index,index)+k;
end
q=1000;
A=b*a;

CB=false(Ndof,1);
CB([7 8 9])=true;
Kr=kG(~CB,~CB);
T=zeros(Ndof,1);
T(CB)=100;
Td=T(~CB);
Q=ones(Ndof,1)*q*A/4;
Q([4 2 6 8])=Q(4)*2;
Q(5)=Q(4)*2;
Qr=Q(~CB);
Tr=Kr\Qr;
T(~CB)=Tr;
bandplot(meinDOF,nodos, T)