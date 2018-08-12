viga=@(E,I,L)(E*I/L^3)*[12 6*L -12 6*L;
    6*L 4*L^2 -6*L 2*L^2;
    -12 -6*L 12 -6*L;
    6*L 2*L^2 -6*L 4*L^2];
Ne=50;
N=Ne+1;
Ndof=2*N;

E=210e9;%Pa
A=0.5e-2;%m^2
I=0.5e-4;%m^4
L=8;%m
h=L/Ne;
w=300e3;%N/m (300kN/m)
k=viga(E,I,h);

elementos=zeros(Ne,4);
for i=1:Ne
    for j=1:4
        elementos(i,j)=(i-1)*2+j;
    end
end

kG=zeros(Ndof,Ndof);

for i=1:Ne
    kG(elementos(i,:),elementos(i,:))=kG(elementos(i,:),elementos(i,:))+k;
end

q=@(x)-w*x/L;

R=zeros(Ndof,1);

for i=1:2:2*Ne
    x=(i-1)*h/2;
    Atri=h/2*(q(x+h)-q(x));
    Acuadr=h*q(x);
    R(i)=R(i)+Acuadr/2+Atri/3;
    R(i+2)=Acuadr/2+Atri*2/3;
end

for i=2:2:2*Ne
    a=sign(q(x+h)-q(x));%Signo de gradiente
    x=(i-2)*h/2;
    if a>=0
        Mtri=h^2/12*(q(x+h)-q(x));
        Mcuadr=h^2*q(x)/24;
        R(i)=R(i)+Mcuadr+.4*Mtri;
        R(i+2)=-Mcuadr-.6*Mtri;
%         R(i)=R(i)+Acuadr/2+Atri/3;
%         R(i+2)=Acuadr/2+Atri*2/3; 
    else
        Mtri=h^2/12*(q(x)-q(x+h));
        Mcuadr=h^2*q(x+h)/24;
        R(i+2)=R(i)-Mcuadr-.4*Mtri;
        R(i)=Mcuadr+.6*Mtri;
    end
    
end

CB=false(Ndof,1);
CB(1)=true;
CB(2)=true;
CB(end)=true;
CB(end-1)=true;
F=R(~CB);
K=kG(~CB,~CB);
U=K\F;
D=zeros(1,Ndof/2);
k=0;
for i=1:2:length(U)
    k=k+1;
    D(k)=U(i);
end
plot(D);