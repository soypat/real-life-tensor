nod=[0;1500;3000];
elenod=[1 2;2 3];
Ne=size(elenod,1);
N=size(nod,1);
Nd=N*1;
b=20;
h=150;
A=b*h;
E=210e3;
a=1.2e-5;
Le=zeros(Ne,1);
kG=zeros(Nd);
for i=1:Ne
ns=nod(elenod(i,1));
ne=nod(elenod(i,2));
lx=ne-ns;
Le(i)=lx;
klocal=A*E/Le(i)*[1 -1;-1 1];
index=[elenod(i,1) elenod(i,2)];
kG(index,index)=kG(index,index)+klocal;
end
%Cargas termicas
Dt=100;%Celsius
Ft=A*E*a*Dt;

%
P=1e3;
CB=false(N,1);
CB([1 3])=true;
Kr=kG(~CB,~CB);
R=zeros(N,1);
R(2)=P-Ft;
F=R(~CB);
U=Kr\F;
D=zeros(N,1);
D(~CB)=U;
Fext=kG*D
H1=