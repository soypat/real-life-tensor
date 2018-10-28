L=.5;
b=1;
close all
Nx=4; %Nodos en los que dividimos el problema
Ny=4;
iterants=[fibonacci(7) fibonacci(8) fibonacci(9) fibonacci(10) fibonacci(11) fibonacci(12)];
err=zeros(length(iterants),1);
k=0; %contador (no es la conductividad)
for N=iterants
k=k+1;
Nx=N;
Ny=N;
dx=L/(Nx-1);
dy=b/(Ny-1);
x=0:dx:L;
y=0:dy:b;

doftot=Nx*Ny;
DOF=reshape(1:doftot,Ny,Nx);
C=sparse(doftot);
Q=zeros(doftot,1);
T=zeros(doftot,1);

nwall=DOF(1,:);
swall=DOF(end,:);
wwall=DOF(:,1);
ewall=DOF(:,end);
Tm=20;
%Tenemos una c
Tnwall=Tm*sin(pi.*x./L);

C(nwall,nwall)=eye(length(nwall));
C(swall,swall)=eye(length(swall));
C(wwall,wwall)=eye(length(wwall));
C(ewall,ewall)=eye(length(ewall));

Q(nwall)=Tnwall;

for j = 2:Nx-1
    for i= 2:Ny-1
        C(DOF(i,j),DOF(i,j))=C(DOF(i,j),DOF(i,j))-2*(1/dx^2+1/dy^2);
        C(DOF(i,j),DOF(i+1,j))=C(DOF(i,j),DOF(i+1,j))+1/dy^2;
        C(DOF(i,j),DOF(i-1,j))=C(DOF(i,j),DOF(i-1,j))+1/dy^2;
        C(DOF(i,j),DOF(i,j+1))=C(DOF(i,j),DOF(i,j+1))+1/dx^2;
        C(DOF(i,j),DOF(i,j-1))=C(DOF(i,j),DOF(i,j-1))+1/dx^2;        
    end
end
T=C\Q;
Tgrid=reshape(T,Ny,Nx);
Texact=solucion_analitica2(x,y);
err(k)=(numel(Tgrid)^-1)*sum(sum((Tgrid-Texact).^2));
end

% bar3(y,Tgrid)
% subplot(3,1,1)
figure(1)
contourf(x,y,Tgrid,'ShowText','on')
title('Solución obtenida (Isotermas)')
% subplot(2,1,2)
% contourf(x,y,Texact)
% title('Solución Exacta')
figure(2)
subplot(1,3,1)
surf1=surf(x,y,Tgrid,'edgecolor','none');
title('Solución Obtenida')
 view(0,90)
 daspect([1 1 1])
subplot(1,3,2)
surf2=surf(x,y,Texact,'edgecolor','none');
 view(0,90)
 daspect([1 1 1])
title('Solución Exacta')
subplot(1,3,3);
loglog(iterants.^2,err)
title('Error cuadratico medio (Escala logarítmica)')
xlabel('Cantidad de volúmenes');
ylabel('Error')
surf2.FaceColor='interp';
%% Conservacion de energía
% Qnwall=zeros(length(nwall),1)
% for i=nwall(2:end-1)
%     T(DOF()

Qewall= zeros(length(ewall),1);
% Tewall = zeros(length(ewall),1);
Tewall = Tgrid(2:end-1,end-1);
Twwall = Tgrid(2:end-1,2);
Tnwall = Tgrid(2,2:end-1);
Tswall = Tgrid(end-1,2:end-1);



Qwest= -dy*Twwall/dx;
Qeast= -dy*Tewall/dx;
Qsouth= -dx*Tswall/dy;
Qnorth= -dx/dy*(Tnwall-Tgrid(1,2:end-1));
sum(Qwest)
sum(Qnorth)
sum(Qeast)
sum(Qsouth)
dQ=sum(Qwest)+sum(Qnorth)+sum(Qeast)+sum(Qsouth);
% Tnwall = Tgrud

%% dQ exacta

Tewalle = Texact(2:end-1,end-1);
Twwalle = Texact(2:end-1,2);
Tnwalle = Texact(2,2:end-1);
Tswalle = Texact(end-1,2:end-1);
Qweste= -dy*Twwalle/dx;
Qeaste= -dy*Tewalle/dx;
Qsouthe= -dx*Tswalle/dy;
Qnorthe= -dx/dy*(Tnwalle-Texact(1,2:end-1));

dQe=sum(Qweste)+sum(Qnorthe)+sum(Qeaste)+sum(Qsouthe);