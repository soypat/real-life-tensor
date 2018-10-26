E=210e3;
I=2e3;
L=2e3;
k1=(E*I/L^3)*[12 6*L -12 6*L;
    6*L 4*L^2 -6*L 2*L^2;
    -12 -6*L 12 -6*L;
    6*L 2*L^2 -6*L 4*L^2];
k2=k1;
elementos=[1,2,3,4;3,4,5,6];
kG=zeros(6);
kG(elementos(1,:),elementos(1,:))=kG(elementos(1,:),elementos(1,:))+k1;
kG(elementos(2,:),elementos(2,:))=kG(elementos(2,:),elementos(2,:))+k2;

