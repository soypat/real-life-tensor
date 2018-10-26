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
elementos=zeros(Ne,4);

for i=1:Ne
    for j=1:4
        elementos(i,j)=(i-1)*2+j;
    end
end
k1=viga(E,I,L);
k2=hingeL(E,I,1);
k3=viga(E,I,1);
