%VIGORRACIÓN
%Pa  m   N
E=210e9;
I=2e-4;
L=2;

hingeL=@(E,I,L) (3*E*I/L^3)*[1 L -1 0;
    L L^2 -L 0;
    -1 -L 1 0;
    0 0 0 0];
viga=@(E,I,L)(E*I/L^3)*[12 6*L -12 6*L;
    6*L 4*L^2 -6*L 2*L^2;
    -12 -6*L 12 -6*L;
    6*L 2*L^2 -6*L 4*L^2];
Ne=3;
N=Ne+1;
ndof=2;
Ndof=N*ndof;
elementos=zeros(Ne,4);

for i=1:Ne
    for j=1:4
        elementos(i,j)=(i-1)*2+j;
    end
end
k1=viga(E,I,L);
k2=hingeL(E,I,1);
k3=viga(E,I,1);

kG=zeros(Ndof);
for j=1:Ne
    k=eval(sprintf('k%0.0f',j));
    kG(elementos(j,:),elementos(j,:))=kG(elementos(j,:),elementos(j,:))+k;
end
R=zeros(Ndof,1);
w=10e3;
q=@(x) 10e3+x*0;

R(5)=-w/2;
R(7)=-w/2;
R(6)=-w/12;
R(8)=w/12;

CB=false(Ndof,1);
CB([1 2])=[true true];
CB([7 8])=[true true];
CB(3)=true;

F=R(~CB);
K=kG(~CB,~CB);
U=K\F