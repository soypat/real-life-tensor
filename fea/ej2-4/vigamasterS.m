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

E=210e9;  %Pa
I=0.5e-4; %m^4
A=0.5e-2; %m^2

nodos=[0 0;
    0 4;
    10 4;
    0 8;
    10 8];
elementos=[1 2 3 4 5 6;
    4 5 6 7 8 9;
    4 5 6 10 11 12;
    10 11 12 13 14 15];
q=@(y) 300e3*(8-y)/8;

Ne=4;
N=5;
ndof=3;
Ndof=ndof*N;

kvert=vigota(E,A,I,4);
khor=vigota(E,A,I,10);
k1=Tvu(90)'*kvert*Tvu(90);
k2=Tvu(0)'*kvert*Tvu(0);
k3=k1;
k4=k2;

kG=zeros(Ndof);

kG(elementos(1,:),elementos(1,:))=kG(elementos(1,:),elementos(1,:))+k1;
kG(elementos(2,:),elementos(2,:))=kG(elementos(2,:),elementos(2,:))+k2;
kG(elementos(3,:),elementos(3,:))=kG(elementos(3,:),elementos(3,:))+k3;
kG(elementos(4,:),elementos(4,:))=kG(elementos(4,:),elementos(4,:))+k4;

R=zeros(Ndof,1);
h=4;
elemento=1;
R=qapplyx(q,0,h,1,2,R,ndof);
R=qapplyx(q,4,h,2,4,R,ndof);

CB=false(Ndof,1);
CB([1 2 3])=[true true true];
CB([7 8 9])=[true true true];
CB([13 14 15])=[true true true];

K=kG(~CB,~CB);
F=R(~CB);
U=K\F;


