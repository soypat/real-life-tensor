%%Dimensiones del problema
clear
convection=false;
L=1; % metros
Area=.1; %m^2 es la sección
k=1; %Watt per meter kelvin
D=sqrt(4*Area/pi); 
qvol=10; %W/m^3

%% Condiciones Borde
T0=0;
Tamb=30;
hamb=10; %W/m^2K
%% Begin problem
Ne=10; %Elementos
N=Ne+1;
Le=L/Ne; %long. elementos
nodos=0:Le:L
C=zeros(N);
Q=zeros(N,1);

%volumen 1
C(1,1)=-3/2*k*Area/Le-hamb*pi*D*Le;
C(1,2)=k*Area/Le;
%Revisar

%volumen N 
C(end,end-1)=k*Area/Le;
C(end,end)=-k*Area/Le-hamb*pi*D*Le; %Caso flujo fijado.

for i=2:Ne-1
    % No consideramos calor almacenado (inercia termica) d/dt=0 R.E.
    C(i,i-1)=k*Area/Le;
    C(i,i)=-2*k*Area/Le-hamb*pi*D*Le;
    C(i,i+1)=k*Area/Le;
end

Q(2:end)=-qvol*Area*Le-hamb*pi*D*Le*Tamb;

T=C\Q;
midpoint=@(X) (X(1:end-1)+X(2:end))/2;
figure(1)
plot(nodos,T,'*',[0 nodos],solucion_analitica1([0,nodos]))
title('Perfil de Temperaturas')
ylabel('Temperatura')
xlabel('Posicion')