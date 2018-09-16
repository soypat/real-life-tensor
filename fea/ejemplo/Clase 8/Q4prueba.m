%% Práctica 8 el culo te abrocho
clear; close all; clc
%% Datos
E = 5E9; %[Pa]
nu = 0.3;
lam = E*nu/(1+nu)/(1-2*nu);
mu = E/2/(1+nu);
t = 1;
%% Nodos y elementos
nod = load('nodos.txt');
nod = nod(:,[2 3]);
nnod = size(nod,1);
elem = load('elementos.txt');
elem = elem(:,2:5);
nelem = size(elem,1);

%% Funciones de forma
N1 =@(x,y) (1-x)*(1-y)/4;
N2 =@(x,y) (1+x)*(1-y)/4;
N3 =@(x,y) (1+x)*(1+y)/4;
N4 =@(x,y) (1-x)*(1+y)/4;

%% DOF
dofpornodo = 2;
nodporelem = 4;
doftot = dofpornodo*nnod;
dof = reshape(1:doftot,dofpornodo,nnod)';

%% Matriz constitutiva
% Plane strain
C = [1-nu  nu  0
     nu 1-nu 0
     0 0 1-2*nu]*lam/nu;
 
%% Puntos de gauss 2x2
a = 1/sqrt(3);
upg = [-a -a
       a -a
       a a
       -a a];
npg = size(upg,1);
wpg = ones(npg,1);

%% Matriz
Kglobal = zeros(doftot);
for e = 1:nelem
    Ke = zeros(dofpornodo*nodporelem);
    nodelem = nod(elem(e,:),:);
    for ig = 1:npg
        xi = upg(ig,1);
        eta = upg(ig,2);
        dNxe = [-(1-eta) 1-eta 1+eta -(1+eta)
                   -(1-xi) -(1+xi) 1+xi 1-xi]/4;
        J = dNxe*nodelem;
        dNxy = J\dNxe;
        B = zeros(size(C,2),dofpornodo*nodporelem);
        B(1,1:2:7) = dNxy(1,:);
        B(2,2:2:8) = dNxy(2,:);
        B(3,1:2:7) = dNxy(2,:);
        B(3,2:2:8) = dNxy(1,:);
        Ke = Ke + B'*C*B*t*wpg(ig)*det(J);
    end
    dofs = reshape(dof(elem(e,:),:)',[],1);
    Kglobal(dofs,dofs) = Kglobal(dofs,dofs) + Ke;

    
end
Kglobal = sparse(Kglobal);

%% Cargas de volumen
P = zeros(doftot/2,2);
f1 = 2*(mu+lam)*10^-2; % [N/m^3]
f2 = -(mu+lam)*10^-2; % [N/m^3]
for e=1:nelem
    nodelem = nod(elem(e,:),:);
    for ig= 1:npg
        xi = upg(ig,1);
        eta = upg(ig,2);
        dNxe = [-(1-eta) 1-eta 1+eta -(1+eta)
                -(1-xi) -(1+xi) 1+xi 1-xi]/4;
        J = dNxe*nodelem;
        N = [N1(xi,eta) N2(xi,eta) N3(xi,eta) N4(xi,eta)];
        P(elem(e,:),1) = P(elem(e,:),1) + N'*f1*t*wpg(ig)*det(J);
        P(elem(e,:),2) = P(elem(e,:),2) + N'*f2*t*wpg(ig)*det(J);
    end
end

%% Puntos de gauss 2
a = sqrt(0.6);
upg = [-a a];
npg = size(upg,1);
wpg = [1 1];

%% Cargas de superficie
g11 =@(y) ((2*mu+lam)*y-2*lam)*10^-2; %[N/m^2]
g12 =@(y) mu*(1-2*y)*10^-2; %[N/m^2]
g21 =@(x) mu*(x-2)*10^-2; %[N/m^2]
g22 =@(x) (-2*(2*mu+lam)*x+lam)*10^-2; %[N/m^2]
for e = 1:nelem
    nodelem = nod(elem(e,:),:);
    if nodelem(2,1)==1 % Significa que estoy parado en la cara derecha
        for ig = 1:npg
            xi = 1;
            eta = upg(ig);
            dNxe = [-(1-eta) 1-eta 1+eta -(1+eta)
                    -(1-xi) -(1+xi) 1+xi 1-xi]/4;
            J = dNxe*nodelem;
            delta = (nodelem(3,2)-nodelem(2,2))/2*(1-a);
            sig = [N2(xi,eta) N3(xi,eta)]*[g11(nodelem(2,2)+delta) g11(nodelem(3,2)-delta)]';
            tau = [N2(xi,eta) N3(xi,eta)]*[g12(nodelem(2,2)+delta) g12(nodelem(3,2)-delta)]';
            P(elem(e,2),1) = P(elem(e,2),1) + N2(xi,eta)*sig*J(1,1)*t*wpg(ig);
            P(elem(e,2),2) = P(elem(e,2),2) + N2(xi,eta)*tau*J(1,1)*t*wpg(ig);
            P(elem(e,3),1) = P(elem(e,3),1) + N3(xi,eta)*sig*J(1,1)*t*wpg(ig);
            P(elem(e,3),2) = P(elem(e,3),2) + N3(xi,eta)*tau*J(1,1)*t*wpg(ig);
        end
    end
    if nodelem(3,2)==1 %significa que estoy parado en la cara de arriba
        for ig = 1:npg
            xi = upg(ig);
            eta = 1;
            dNxe = [-(1-eta) 1-eta 1+eta -(1+eta)
                    -(1-xi) -(1+xi) 1+xi 1-xi]/4;
            J = dNxe*nodelem;
            delta = (nodelem(3,1)-nodelem(4,1))/2*(1-a);
            sig = [N3(xi,eta) N4(xi,eta)]*[g22(nodelem(3,1)-delta) g22(nodelem(4,1)+delta)]';
            tau = [N3(xi,eta) N4(xi,eta)]*[g21(nodelem(3,1)-delta) g21(nodelem(4,1)+delta)]';
            P(elem(e,4),1) = P(elem(e,4),1) + N4(xi,eta)*tau*J(1,1)*t*wpg(ig);
            P(elem(e,4),2) = P(elem(e,4),2) + N4(xi,eta)*sig*J(1,1)*t*wpg(ig);
            P(elem(e,3),1) = P(elem(e,3),1) + N3(xi,eta)*tau*J(1,1)*t*wpg(ig);
            P(elem(e,3),2) = P(elem(e,3),2) + N3(xi,eta)*sig*J(1,1)*t*wpg(ig);
        end
    end
end
P = reshape(P',[],1);

%% BC
fijo = zeros(doftot/2,2);
fijo(nod(:,2)==0,:) = true;
fijo(nod(:,1)==0,:) = true;
fijo = reshape(logical(fijo)',[],1);
libre = ~fijo;

%% Reducción
Kglobalred = Kglobal(libre,libre);
Pred = P(libre);

%% Solver
Dred = Kglobalred\Pred;
D = zeros(doftot,1);
D(libre) = Dred;
meshplot(elem,nod,'r')
hold on
nodfinal = nod+reshape(D,dofpornodo,nnod)';
meshplot(elem,nodfinal,'b')
















