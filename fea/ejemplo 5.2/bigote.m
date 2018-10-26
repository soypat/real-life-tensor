vigota=@(E,A,I,L) [A*E/L 0 0 -A*E/L 0 0;
    0 12*E*I/L^3 6*E*I/L^2 0 -12*E*I/L^3 6*E*I/L^2;
    0 6*E*I/L^2 4*E*I/L 0 -6*E*I/L^2 2*E*I/L;
    -A*E/L 0 0 A*E/L 0 0;
    0 -12*E*I/L^3 -6*E*I/L^2 0 12*E*I/L^3 -6*E*I/L^2;
    0 6*E*I/L^2 2*E*I/L 0 -6*E*I/L^2 4*E*I/L];
Tvu=@(phi) [cosd(phi) sind(phi) 0 0 0 0;
    -sind(phi) cosd(phi) 0 0 0 0;
    0 0 1 0 0 0;
    0 0 0 cosd(phi) sind(phi) 0;
    0 0 0 -sind(phi) cosd(phi) 0;
    0 0 0 0 0 1];

elementos=[1 2 3 4 5 6;
    4 5 6 7 8 9];
E=30e6;%psi
A=100;%in
I=1e3;%in
q=@(x)1000/12;
L1=sqrt(30^2+30^2)*12;%in
L2=40*12;
k1=Tvu(45)'*vigota(E,A,I,L1)*Tvu(45);
k2=Tvu(0)'*vigota(E,A,I,L2)*Tvu(0);
ndof=3;
N=3;
Ne=2;
Ndof=ndof*N;


kG=zeros(Ndof);

kG(elementos(1,:),elementos(1,:))=kG(elementos(1,:),elementos(1,:))+k1;
kG(elementos(2,:),elementos(2,:))=kG(elementos(2,:),elementos(2,:))+k2;

F=[0; -20e3; -1000/12*L2^2/12];

CB=false(Ndof,1);
CB([1 2 3 7 8 9])=[true true true true true true];

K=kG(~CB,~CB);
U=K\F