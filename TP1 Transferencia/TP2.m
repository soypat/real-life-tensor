k=1;
L=0.5;
b=1;

Nx=30; %Nodos en los que dividimos el problema
Ny=30;

dx=L/(Nx-1);
dy=b/(Ny-1);
x=0:dx:L;
y=0:dy:b;

doftot=Nx*Ny;
DOF=reshape(1:doftot,Ny,Nx);
C=zeros(doftot);
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
        C(DOF(i,j),DOF(i,j))=C(DOF(i,j),DOF(i,j))-2/dx^2-2/dy^2;
        C(DOF(i,j),DOF(i+1,j))=C(DOF(i,j),DOF(i+1,j))+1/dx^2;
        C(DOF(i,j),DOF(i-1,j))=C(DOF(i,j),DOF(i-1,j))+1/dx^2;
        C(DOF(i,j),DOF(i,j+1))=C(DOF(i,j),DOF(i,j+1))+1/dy^2;
        C(DOF(i,j),DOF(i,j-1))=C(DOF(i,j),DOF(i,j-1))+1/dy^2;        
    end
end

T=C\Q;
Tgrid=reshape(T,Ny,Nx);

% bar3(y,Tgrid)
surf(x,y,Tgrid)










