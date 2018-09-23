%% Elementos Q4
syms x y real
b=1; %ancho dividido 2
h=.5; %Alto dividido 2
a=h;
N1s=(b-x)*(h-y)/(4*b*h); %Eje centrado
N2s=(b+x)*(h-y)/(4*b*h);
N3s= (b+x)*(h+y)/(4*b*h);
N4s= (b-x)*(h+y)/(4*b*h);%Lo saque del Daryl Logan. Está tambien en el Cook
N=[N1s N2s N3s N4s]; %Una fila de matriz de funciones de forma asi puedo derivar split
B=[diff(N,x);diff(N,y)];
C=[25 0;0 25]; %25 W/mk isotropo
k=int(int(B'*C*B,x,-a,a),y,-b,b )
eig(k) % esta es lo que debería dar! Si!
%% Comenzamos a buscar la solucion a problema elementos finitos
meinDOF=[1 2 4 5;2 3 5 6;4 5 7 8;5 6 8 9];
Ne=4
Ndof=9;
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
Kr\Qr