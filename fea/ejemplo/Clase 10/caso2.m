%% Carga distribuida en elementos isoparamétricos
clear; close all; clc
%% Datos
E = 210E9; %[Pa]
nu = 0.3;
lam = E*nu/(1+nu)/(1-2*nu);
mu = E/2/(1+nu);
t = 1;
eleT = 'Q4';

%% Nodos y elementos
nod = [0 0
       8 0
       3 4
       8 4]; % [m]
nnod = size(nod,1);
elem = [1 2 4 3];
nelem = size(elem,1);
nodporelem = 4;

%% Funciones de forma
N1 =@(x,y) (1-x)*(1-y)/4;
N2 =@(x,y) (1+x)*(1-y)/4;
N3 =@(x,y) (1+x)*(1+y)/4;
N4 =@(x,y) (1-x)*(1+y)/4;

%% DOF
dofpornodo = 2;
doftot = dofpornodo*nnod;
dof = reshape((1:doftot)',dofpornodo,nnod)';

%% Constitutiva
% Plane stress
C = [1 nu 0
     nu 1 0
     0 0 1-nu]*E/(1-nu^2);
 
%% Matriz
[wpg, upg, npg] = gauss([2 2]); 
Kglobal = zeros(doftot);
for e = 1:nelem
    nodelem = nod(elem(e,:),:);
    Ke = zeros(nodporelem*dofpornodo);
    for ipg = 1:npg
        ksi = upg(ipg,1);
        eta = upg(ipg,2);
        dNke = shapefunsder([ksi eta],eleT);
        J = dNke*nodelem;
        dNxy = J\dNke;
        B = zeros(size(C,2),nodporelem*dofpornodo);
        B(1,1:2:7) = dNxy(1,:);
        B(2,2:2:8) = dNxy(2,:);
        B(3,1:2:7) = dNxy(2,:);
        B(3,2:2:8) = dNxy(1,:);
        Ke = Ke + B'*C*B*t*wpg(ipg)*det(J);
    end
    dofs = reshape(dof(elem(e,:),:)',[],1);
    Kglobal(dofs,dofs) = Kglobal(dofs,dofs) + Ke;
end

%% Cargas
P = zeros(doftot/2,2);
f =@(x) 5000000000*(1-x/3);
tita = atan(3/4);
npg = 2;
[wpg, upg] = gauss1D(npg);
for e = 1:nelem
    nodelem = nod(elem(e,:),:);
    if nodelem(4,1)/nodelem(4,2) == 3/4
%         for ig = 1:npg
%             ksi = -1;
%             eta = upg(ig);
%             dNke = shapefunsder([ksi eta],eleT);
%             J = dNke*nodelem;
%             delta = (nodelem(4,1)-nodelem(1,1))/2;
%             centronodo = (nodelem(4,1)+nodelem(1,1))/2;
%             phiy = f(centronodo+delta*eta);
%             sig = phiy*sin(tita);
%             tau = phiy*cos(tita);
%             P(elem(e,1),:) = P(elem(e,1),:) + N1(ksi,eta)*wpg(ig)*t*[tau*J(2,1)-sig*J(2,2) tau*J(2,2)+sig*J(2,1)];
%             P(elem(e,4),:) = P(elem(e,4),:) + N4(ksi,eta)*wpg(ig)*t*[tau*J(2,1)-sig*J(2,2) tau*J(2,2)+sig*J(2,1)];
%         end
        for ig = 1:npg
            ksi = -1;
            eta = upg(ig);
            dNke = shapefunsder([ksi eta],eleT);
            J = dNke*nodelem;
            phiy = [f(nodelem(1,1)) f(nodelem(4,1))]';
            N = [N1(ksi,eta) N4(ksi,eta)];
            P(elem(e,[1 4]),2) = P(elem(e,[1 4]),2) + wpg(ig)*t*N'*N*phiy*J(2,1);
        end
    end
end
P = reshape(P',[],1);

%% BC
fijo = zeros(doftot/2,2);
fijo(nod(:,1)==8,:) = true;
fijo = reshape(fijo',[],1);
libre = ~fijo;

%% Solver
Dred = Kglobal(libre,libre)\P(libre);
D = zeros(doftot,1);
D(libre) = Dred;
meshplot(elem,nod,'r')
hold on
nodfinal = nod+reshape(D,dofpornodo,nnod)';
meshplot(elem,nodfinal,'b')

%% Tensiones
stress = zeros(nelem,nodporelem,3);
unod = [ -1 -1
          1 -1
          1  1
         -1  1 ];
for e = 1:nelem
    nodelem = nod(elem(e,:),:);
    for in = 1:nodporelem
        % Punto de Gauss
        ksi = unod(in,1);
        eta = unod(in,2);  
        % Derivadas de las funciones de forma respecto de ksi, eta
        dNke = shapefunsder([ksi eta],eleT);
        % Derivadas de x,y, respecto de ksi, eta
        J = dNke*nodelem;                      
        % Derivadas de las funciones de forma respecto de x,y.
        dNxy = J\dNke;
        B = zeros(size(C,2),dofpornodo*nodporelem);
        B(1,1:2:7) = dNxy(1,:);
        B(2,2:2:8) = dNxy(2,:);
        B(3,1:2:7) = dNxy(2,:);
        B(3,2:2:8) = dNxy(1,:);
        dofs = reshape(dof(elem(e,:),:)',[],1);
        stress(e,in,:) = C*B*D(dofs);
    end
end

bandplot(elem,nodfinal,stress(:,:,1),[],'k')












