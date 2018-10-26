E=210e3; %MPa=N/mm^2
A=5e4;%mm^2
k=2e3; %N/mm Simetria, la aprovecho
thetad=60;%de=90
nine=0;
P=50e3;%N simetria
L=5e3;%mm

N=2;

Ne=N-1;
h=L/Ne;
k1=sind(thetad)^2*E*A/L;
% nodos=[0:h:L];
% elementos=zeros(Ne,2);
% for i=1:Ne
%     elementos(i,1)=i;
%     elementos(i,2)=i+1;
% end
elementos=[1 2 3 4;
            5 6 7 8;
            9 10 11 12;
            13 14 15 16];

T1=[cosd(nine) 0;
    sind(nine) 0;
    0 cosd(nine);
    0 sind(nine)];
T2=[cosd(thetad) 0;
    sind(thetad) 0;
    0 cosd(thetad);
    0 sind(thetad)];
km1=[k -k;
    -k k];
km2=[k1 -k1;
    -k1 k1];
km1=T1*km1*T1';
km2=T2*km2*T2';
kG2=zeros(6);
kG2(elementos(1,:),elementos(1,:))=kG2(elementos(1,:),elementos(1,:))+km1;
kG2(elementos(3,:),elementos(3,:))=kG2(elementos(3,:),elementos(3,:))+km2;
kG=[k -k 0;
    -k k+k1 -k1; 
    0 -k1 k1];

K=kG(2:end-1,2:end-1); %Condiciones de borde

f=[0;
    -P;
    0];
R=f(2:end-1);
U=K\R;
