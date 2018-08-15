nodos=12*[-5,0;-5,10;5,10;5,0];%in
E=30e6;%psi
A=10;%in^2
I13=200;%in^4
I2=100;
L=10*12;
% viga=@(E,I,L)(E*I/L^3)*[12 6*L -12 6*L;
%     6*L 4*L^2 -6*L 2*L^2;
%     -12 -6*L 12 -6*L;
%     6*L 2*L^2 -6*L 4*L^2];
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
Ne=3;
N=Ne+1;
ndof=3; %dof unitarios
Ndof=N*ndof;
elementos=zeros(Ne,6);
for i = 1:Ne
    for j = 1:6
        elementos(i,j)=(i-1)*ndof+j;
    end
end



kG=zeros(Ndof);

T90=Tvu(90);
T0=Tvu(0);
Tm90=Tvu(-90);

k13=vigota(E,A,I13,L);
k=vigota(E,A,I2,L);
k1=T90'*k13*T90;
k2=T0'*k*T0;
k3=Tm90'*k13*Tm90;

kG(elementos(1,:),elementos(1,:))=kG(elementos(1,:),elementos(1,:))+k1;
kG(elementos(2,:),elementos(2,:))=kG(elementos(2,:),elementos(2,:))+k2;
kG(elementos(3,:),elementos(3,:))=kG(elementos(3,:),elementos(3,:))+k3;

CB=false(Ndof,1);
CB(1)=true;
CB(2)=true;
CB(3)=true;
CB(end-2)=true;
CB(end-1)=true;
CB(end)=true;

K=kG(~CB,~CB);

R=zeros(Ndof,1);
R(4)=10e3;
R(9)=5e3;
F=R(~CB);
U=K\F