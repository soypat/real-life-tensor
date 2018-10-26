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
Nx=6; %cantidad de elementos (2 medios volumenes y Nx-2 completos )
dx=L/(Nx-1); %long. elementos
nodos=0:dx:L;
C=zeros(Nx);
Q=zeros(Nx,1);

%volumen 1
C(1,1)=1;
Q(1)=T0;
%volumen N 
C(end,end-1)=1/2*k*Area/dx;
C(end,end)=1/2*(-k*Area/dx-hamb*pi*D*dx); %Caso flujo fijado.
Q(end)=1/2*(-qvol*Area*dx-hamb*pi*D*dx*Tamb);

for i=2:Nx-1
    % No consideramos calor almacenado (inercia termica) d/dt=0 R.E.
    C(i,i-1)=k*Area/dx;
    C(i,i)=-2*k*Area/dx-hamb*pi*D*dx;
    C(i,i+1)=k*Area/dx;
    Q(i)=-qvol*Area*dx-hamb*pi*D*dx*Tamb;
end

T=C\Q;
figure(1)
hold on
x_sol=0:0.05:1;
plot(x_sol,solucion_analitica1(x_sol),'b')
plot(nodos,T,'k-.*')
legend('analitic','finite volume')
hold off

title('Perfil de Temperaturas')
ylabel('Temperatura')
xlabel('Posicion')