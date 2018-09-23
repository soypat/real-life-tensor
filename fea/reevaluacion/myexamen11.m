b=0.01;
L=1;
E=79e9;%Pa
A=b^2;
Iz=b*b^3/12;
kele=Kv(E,A,Iz,L);

D=@(u,v,g) [0;v;g;u;-v;-g];
%% VIGA (fixed u1) theta1  |---------------------- u2 theta2
%Let theta1=-theta2
giro=8.45*9e-7; %en radianes
u=1e-4;% en metros
%% Nodo 1 contiene desplazamientos solo en theta, si no no hay fuerza axil y no puedo obtener resultados diferentes.
Dab=D(u,0,giro); %Hay desplazamientos en u y theta en el nodo 2,
x=L/2;
Mz=E*Iz*( (-6/L^2 + 12*x/L^3)*Dab(2)+(-4/L+6*x/L^2)*Dab(3)+(6/L^2-12*x/L^3)*Dab(5)+(-2/L+6*x/L^2)*Dab(6) )
Nx=Dab(1)*-1/L+Dab(4)*1/L
sigma=Nx/A+(Mz*(15/2)/Iz)/1e6 %En MPa un extremo 
sigma=Nx/A-(Mz*(15/2)/Iz)/1e6 %Otro extremo