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
Nx=10; %cantidad de elementos (2 medios volumenes y Nx-2 completos )
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

Texact=solucion_analitica1(0:dx:1)';
error_max=max(abs(T(2:end)-Texact(2:end))./Texact(2:end))*100;
figure(1)
subplot(2,1,1);
plot(x_sol,solucion_analitica1(x_sol),'b')
hold on
plot(nodos,T,'k-.*')
legend('analitic','finite volume')
text(.2,15,sprintf('Error Relativo = %f%%',error_max))
hold off

title('Perfil de Temperaturas')
ylabel('Temperatura')
xlabel('Posicion')

%% Analisis de Conservation of Energy
Qgen=qvol*Area*dx*ones(Nx,1);
Qgen(1)=Qgen(1)/2;
Qgen(Nx)=Qgen(Nx)/2;


Qeast=k*Area*(T(end)-T(end-1))/dx;
Qwest=-k*Area*(T(2)-T(1))/(dx);

Qconv=zeros(Nx,1);
Qconv(1)=hamb*D*pi*dx*(Tamb-T(1))/2;
Qconv(Nx)=hamb*D*pi*dx*(Tamb-T(Nx))/2;
for i=2:Nx-1
    Qconv(i)=hamb*D*pi*dx*(Tamb-T(i));
end

subplot(2,1,2);
plot(nodos,Qconv)
title('Calor convectado sobre barra')
ylabel('Calor Convectado [W/m]')
xlabel('Posicion [m]')
dQ=Qeast+sum(Qconv)+sum(Qgen)+Qwest;
% texty=sprintf('Calor almacenado = %0.3f',dQ);
text(.22,max(Qconv)/2,sprintf('Calor almacenado = %0.3f W',dQ))