%% Caso 3, Q8 con corte tipo piramide de guiza, no de guiso
clear; close all; clc
%% Datos
E = 210E9; %[Pa]
nu = 0.3;
lam = E*nu/(1+nu)/(1-2*nu);
mu = E/2/(1+nu);
t = 1;
eleT = 'Q8';
cargaT = 'Superficie';
%% Nodos y elementos
d = .2; %[m]
nod = [0 0
       1 0
       2 0
       0 1
       2 1
       0 2
       1 2+d
       2 2]; % [m]
nnod = size(nod,1);
elem = [1 3 8 6 2 5 7 4];
nelem = size(elem,1);
nodporelem = 8;

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
[wpg, upg, npg] = gauss([3 3]); 
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
        B(1,1:2:15) = dNxy(1,:);
        B(2,2:2:16) = dNxy(2,:);
        B(3,1:2:15) = dNxy(2,:);
        B(3,2:2:16) = dNxy(1,:);
        Ke = Ke + B'*C*B*t*wpg(ipg)*det(J);
    end
    dofs = reshape(dof(elem(e,:),:)',[],1);
    Kglobal(dofs,dofs) = Kglobal(dofs,dofs) + Ke;
end

%% Cargas
switch cargaT
%% Cargas
    case 'Superficie'
P = zeros(doftot/2,2);
T = 10000000000; %[Pa]
% tauizq =@(x) T*(x+1)/2;
% tauder =@(x) T*(1-x)/2;
tauizq =@(x) T*x;
tauder =@(x) T*(2-x);
npg = 2;
[wpg, upg] = gauss1D(npg);
for e = 1:nelem
    nodelem = nod(elem(e,:),:);
%     if nodelem(4,:) == [0, 2]
%         for ig = 1:npg %izquierda
%             ksi = upg(ig);
%             eta = 1;
%             N = shapefuns([ksi,eta],eleT);
%             dNke = shapefunsder([ksi eta],eleT);
%             J = dNke*nodelem;
%             tau = tauizq(ksi);
%             P(elem(e,4),:) = P(elem(e,4),:) + N(4)*wpg(ig)*t*[tau*J(1,1) tau*J(1,2)];
%             P(elem(e,7),:) = P(elem(e,7),:) + N(7)*wpg(ig)*t*[tau*J(1,1) tau*J(1,2)];
%         end
%         for ig = 1:npg %derecha
%             ksi = upg(ig);
%             eta = 1;
%             N = shapefuns([ksi,eta],eleT);
%             dNke = shapefunsder([ksi eta],eleT);
%             J = dNke*nodelem;
%             tau = tauder(ksi);
%             P(elem(e,7),:) = P(elem(e,7),:) + N(7)*wpg(ig)*t*[tau*J(1,1) tau*J(1,2)];
%             P(elem(e,3),:) = P(elem(e,3),:) + N(3)*wpg(ig)*t*[tau*J(1,1) tau*J(1,2)];
%         end
%     end
    if nodelem(4,:) == [0, 2]
        for ig = 1:npg 
            ksi = upg(ig);
            eta = 1;
            N = shapefuns([ksi eta],eleT);
            dNke = shapefunsder([ksi eta],eleT);
            J = dNke*nodelem;
            tau = [N(4) N(7) N(3)]*[tauizq(nodelem(4,1)) tauizq(nodelem(7,1)) tauder(nodelem(3,1))]';
            P(elem(e,4),:) = P(elem(e,4),:) + N(4)*wpg(ig)*t*[tau*J(1,1) tau*J(1,2)];
            P(elem(e,7),:) = P(elem(e,7),:) + N(7)*wpg(ig)*t*[tau*J(1,1) tau*J(1,2)];
            P(elem(e,3),:) = P(elem(e,3),:) + N(3)*wpg(ig)*t*[tau*J(1,1) tau*J(1,2)];
        end
    end
end
P = reshape(P',[],1);


%% Cargas volumétricas
    case 'Volumnen'
P = zeros(doftot/2,2);
f = 50000; % [N/m^3]
[wpg, upg, npg] = gauss([3 3]); 
for e=1:nelem
    nodelem = nod(elem(e,:),:);
    for ig= 1:npg
        ksi = upg(ig,1);
        eta = upg(ig,2);
        dNxe = shapefunsder([ksi eta],eleT)
        J = dNxe*nodelem;
        N = [N1(ksi,eta) N2(ksi,eta) N3(ksi,eta) N4(ksi,eta)];
        P(elem(e,:),1) = P(elem(e,:),1) + N'*f1*t*wpg(ig)*det(J);
        P(elem(e,:),2) = P(elem(e,:),2) + N'*f2*t*wpg(ig)*det(J);
    end
end
P = reshape(P',[],1);
end

%% BC
fijo = zeros(doftot/2,2);
fijo([1 3],:) = true;
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
         -1  1
          0 -1
          1  0
          0  1
         -1  0];
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
        B(1,1:2:15) = dNxy(1,:);
        B(2,2:2:16) = dNxy(2,:);
        B(3,1:2:15) = dNxy(2,:);
        B(3,2:2:16) = dNxy(1,:);
        dofs = reshape(dof(elem(e,:),:)',[],1);
        stress(e,in,:) = C*B*D(dofs);
    end
end

bandplot(elem,nodfinal,stress(:,:,1),[],'k')












