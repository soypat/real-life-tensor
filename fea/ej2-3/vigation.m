Ne=2;
E=210e9;%Pa
L1=12;%m
L2=6;%m

viga=@(E,I,L)(E*I/L^3)*[12 6*L -12 6*L;
    6*L 4*L^2 -6*L 2*L^2;
    -12 -6*L 12 -6*L;
    6*L 2*L^2 -6*L 4*L^2];


Im=700e-6;%m^4
I=Im; %*(1000)^4; %m^4


k2=viga(E,2*I,L2);
h=L1/Ne;
k1=viga(E,I,h);
nodos1=[0:h:L1];
nodos=[nodos1 L1+L2];

N=Ne+1;
Ndof=4+2*(Ne-1);%dof sobre carga distribuida
Ndoft=Ndof+2;%dof GLOBALES
elementos=zeros(Ne,4);
for i=1:Ne+1
    for j=1:4
        elementos(i,j)=1+(i-1)*2+j-1;
    end
end

kG=zeros(Ndoft,Ndoft);

for i=1:Ne
    kG(elementos(i,:),elementos(i,:))=kG(elementos(i,:),elementos(i,:))+k1;
end
i=i+1;
kG(elementos(i,:),elementos(i,:))=kG(elementos(i,:),elementos(i,:))+k2;
q=@(x)-10e3;
R=zeros(Ndoft,1);

for i=1:2:2*Ne
    x=(i-1)*h;
    Atri=h/2*(q(x+h)-q(x));
    Acuadr=h*q(x);
    R(i)=R(i)+Acuadr/2+Atri/3;
    R(i+2)=Acuadr/2+Atri*2/3;
end
R(end-1)=-100e3;
CB=false(Ndoft,1);
CB(1)=1;%Deflexion```q1 
CB(2*N-1)=1;%Deflexion
CB(2)=1;%Angulo
F=R(~CB);
F(1)=-q(rand())*L1^2/12;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K=kG(~CB,~CB);
U=K\F;
Ugraf=[0; 0; U(1:2*N-4); 0; U(2*N-3:end)];
% plotnodes=zeros(1,N);
D=zeros(Ndoft/2,1);
theta=zeros(Ndoft/2,1);
for i=1:floor(Ndoft/2)
    D(i)=Ugraf(i*2-1);
    theta(i)=Ugraf(i*2);
end


