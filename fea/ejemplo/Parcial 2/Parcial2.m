%% Parcial 2
clear; close all; clc
%% Ej1.a
syms x y real
X = [1 x y x^2 x*y];
A = [1 -1 -1  1  1
     1  1 -1  1 -1
     1  1  1  1  1
     1 -1  1  1 -1
     1  0  0  0  0];
N = X/A;
dNx = diff(N,x);
dNy = diff(N,y);
dN = [dNx; dNy];
sum(N);  % == 1
sum(dN); % == 0
 
%% Datos
E = 200E3; %[MPa]
nu = 0.3;
lam = E*nu/(1+nu)/(1-2*nu);
mu = E/2/(1+nu);
t = 1;

%% Nodos y Elementos
nod = [-1 -1
        1 -1
        0  0
       -1  1
        1  1]*1E3;
nnod = size(nod,1);
elem = [1 2 5 4 3];
nelem = size(elem,1);

%% DOF
dofpornodo = 2;
nodporelem = 5;
doftot = dofpornodo*nnod;
dof = reshape((1:doftot)',dofpornodo,nnod)';

%% Constitutiva
Cstrain = E/((1 + nu)*(1 - 2*nu))*[ 1 - nu      nu            0.0
                                      nu       1 - nu         0.0
                                      0.0        0.0     (1 - 2*nu)/2 ];

%% Gauss para regla de 2x2 (numeración de nodos como figura 6.3-3 pág 212)
a   = 1/sqrt(3);
% Ubicaciones puntos de Gauss
upg = [ -a  -a
        -a   a
         a  -a
         a   a ];    
% Número de puntos de Gauss
npg = size(upg,1);
wpg = ones(npg,1);

%% Matriz de rigidez
Kglobal = zeros(doftot);
for e = 1:nelem
    Ke = zeros(dofpornodo*nodporelem);
    nodelem = nod(elem(e,:),:);
    for ipg = 1:npg
        % Punto de Gauss
        x = upg(ipg,1);
        y = upg(ipg,2);  
        % Derivadas de las funciones de forma respecto de ksi, eta
        dN = [ x/2 + y/4 - 1/4, x/2 - y/4 + 1/4, x/2 + y/4 + 1/4, x/2 - y/4 - 1/4, -2*x
                     x/4 - 1/4,     - x/4 - 1/4,       x/4 + 1/4,       1/4 - x/4,    0];
        % Derivadas de x,y, respecto de ksi, eta
        J = dN*nodelem;                      
        % Derivadas de las funciones de forma respecto de x,y.
        dNxy = J\dN;
        
        B = zeros(size(Cstrain,2),size(Ke,1));
        B(1,1:2:9) = dNxy(1,:);
        B(2,2:2:10) = dNxy(2,:);
        B(3,1:2:9) = dNxy(2,:);
        B(3,2:2:10) = dNxy(1,:); 

        Ke = Ke + B'*Cstrain*B*wpg(ipg)*det(J);
    end
    dofs = dof(elem(e,:),:);
    dofs = reshape(dofs',[],1);
    Kglobal(dofs,dofs) = Kglobal(dofs,dofs) + Ke;  
end

%% BC
fijo = zeros(doftot/2,2);
fijo([1 2],:) = [1 1; 0 1];
fijo = logical(reshape(fijo',[],1));
libre = ~fijo;

%% Cargas
P = zeros(doftot/2,2);
Q = 0.002; %[N/mm]
f =@(x) Q*x;
a   = 1/sqrt(3);
% a = sqrt(0.6);
upg = [-a a];    
npg = size(upg,2);
wpg = ones(npg,1);
% wpg = [5 8 5]/9;
for e = 1:nelem
    nodelem = nod(elem(e,:),:);
    for ipg = 1:npg
        % Punto de Gauss
        x = upg(ipg);
        y = 1;
        N = [ (x*y)/4 - y/4 - x/4 + x^2/4, x/4 - y/4 - (x*y)/4 + x^2/4, x/4 + y/4 + (x*y)/4 + x^2/4, y/4 - x/4 - (x*y)/4 + x^2/4, 1 - x^2];
        dN = [ x/2 + y/4 - 1/4, x/2 - y/4 + 1/4, x/2 + y/4 + 1/4, x/2 - y/4 - 1/4, -2*x
                     x/4 - 1/4,     - x/4 - 1/4,       x/4 + 1/4,       1/4 - x/4,    0];
        sig = f(x);
        J = dN*nodelem;
        P(elem(e,4),:) = P(elem(e,4),:) + N(4)*wpg(ipg)*t*[0 sig*J(1,1)];
        P(elem(e,3),:) = P(elem(e,3),:) + N(3)*wpg(ipg)*t*[0 sig*J(1,1)];
    end
end
P = reshape(P',[],1);

%% Solver
Dred = Kglobal(libre,libre)\P(libre);
D = zeros(doftot,1);
D(libre) = Dred;
cnodfinal = nod+reshape(D,dofpornodo,nnod)'*10000000;
figure
hold on; grid on; axis equal;
plot(nod(:,1),nod(:,2),'.')
plot(cnodfinal(:,1),cnodfinal(:,2),'.')

%% Gauss para regla de 2x2 (numeración de nodos como figura 6.3-3 pág 212)
a   = 1/sqrt(3);
% Ubicaciones puntos de Gauss
upg = [ -a  -a
        -a   a
         a  -a
         a   a ];    
% Número de puntos de Gauss
npg = size(upg,1);
wpg = ones(npg,1);

%% Tensiones
stress = zeros(nelem,nodporelem-1,3);
sigvm = zeros(nelem,nodporelem-1);
for e = 1:nelem
    nodelem = nod(elem(e,1:4),:);
    for in = 1:npg
        x = upg(in,1);
        y = upg(in,2);  
        dN = 1/4*[-(1-y)   1-y    1+y  -(1+y)
                  -(1-x) -(1+x)   1+x    1-x ];
        J = dN*nodelem;                      
        dNxy = J\dN;
        B = zeros(size(Cstrain,2),dofpornodo*npg);
        B(1,1:2:7) = dNxy(1,:);
        B(2,2:2:8) = dNxy(2,:);
        B(3,1:2:7) = dNxy(2,:);
        B(3,2:2:8) = dNxy(1,:);
        dofs = reshape(dof(elem(e,1:4),:)',[],1);
        stress(e,in,:) = Cstrain*B*D(dofs);
        sigmavm(e,in) = sqrt(stress(e,in,1)^2-stress(e,in,1)*stress(e,in,2)+stress(e,in,2)^2+3*stress(e,in,3)^2);
    end
end











