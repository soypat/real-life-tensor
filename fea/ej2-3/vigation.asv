viga=@(E,I,L)(E*I/L^3)*[12 6*L -12 6*L;
    6*L 4*L^2 -6*L 2*L^2;
    -12 -6*L 12 -6*L;
    6*L 2*L^2 -6*L 4*L^2];
E=210e3;
L1=12e3;
L2=6e3;

Im=700e-6;%m^4
I=Im*(1000)^4;

k1=viga(210,I,L1);
k2=viga(210,2*I,L2);

Ne=2;
N=Ne+1;
Ndof=4+2*(Ne-1);
elementos=zeros(Ne,4);
for i=1:Ne+1
    for j=1:4
        elementos(i,j)=1+(i-1)*2+j-1;
    end
end

kG=zeros(Ndof+2);

for i=1:Ne
    kG(elementos(i,:),elementos(i,:))=kG(elementos(i,:),elementos(i,:))+k1;
end
i=i+1;
kG(elementos(i,:),elementos(i,:))=kG(elementos(i,:),elementos(i,:))+k2;
q=@(x)10;
R=zeros(Ndof,1);
for i=1:N
    R(i)=q(1)
end