%% Ej 2
clear; close all; clc

%% Datos
E = 30e9; %[Pa]
NU = 0.18;
rho = 2000; %[kg/m^3]
g = -9.81; %[m/s^2]
t = 1;
%% PROGRAMA DE ELEMENTOS FINITOS PARA ISOPARAMÉTRICOS Q4

elements = load('ElementosT1.txt');

nodes = load('NodosT1.txt');

nDofNod = 2;                    % Número de grados de libertad por nodo
nNodEle = 4;                    % Número de nodos por elemento
nel = size(elements,1);         % Número de elementos
nNod = size(nodes,1);           % Número de nodos
nDofTot = nDofNod*nNod;         % Número de grados de libertad

bc = false(nNod,nDofNod);       % Matriz de condiciones de borde

R = zeros(nNod,nDofNod);        % Vector de cargas

% Propiedades del material


meshplot(elements,nodes,'b')

%% Matriz Constitutiva (plane strain)
 
C = E/((1 + NU)*(1 - 2*NU))*[ 1 - NU      NU            0.0
                                NU       1 - NU         0.0
                               0.0        0.0     (1 - 2*NU)/2 ];
                           
%% Gauss           
a   = 1/sqrt(3);
% Ubicaciones puntos de Gauss
upg = [ -a  -a
         a  -a
         a   a
        -a   a ];    
% Número de puntos de Gauss
npg = size(upg,1);
wpg = ones(npg,1);

%% Matriz de rigidez (Para elementos Q4 con regla de Gauss de 2x2)
K = zeros(nDofTot);
nodeDofs = reshape(1:nDofTot,nDofNod,nNod)';
for iele = 1:nel
    Ke = zeros(nDofNod*nNodEle);
    nodesEle = nodes(elements(iele,:),:);
    for ipg = 1:npg
        % Punto de Gauss
        ksi = upg(ipg,1);
        eta = upg(ipg,2);  
        % Derivadas de las funciones de forma respecto de ksi, eta
        dN = 1/4*[-(1-eta)   1-eta    1+eta  -(1+eta)
                  -(1-ksi) -(1+ksi)   1+ksi    1-ksi ];  
        % Derivadas de x,y, respecto de ksi, eta
        jac = dN*nodesEle;                      
        % Derivadas de las funciones de forma respecto de x,y.
        dNxy = jac\dN;          % dNxy = inv(jac)*dN
        
        B = zeros(size(C,2),nDofNod*nNodEle);
        B(1,1:2:7) = dNxy(1,:);
        B(2,2:2:8) = dNxy(2,:);
        B(3,1:2:7) = dNxy(2,:);
        B(3,2:2:8) = dNxy(1,:); 

        Ke = Ke + B'*C*B*wpg(ipg)*det(jac);
    end
    eleDofs = nodeDofs(elements(iele,:),:);
    eleDofs = reshape(eleDofs',[],1);
    K(eleDofs,eleDofs) = K(eleDofs,eleDofs) + Ke;  
end

%% BC
bc = zeros(nNod,nDofNod);
bc(nodes(:,2)==0,:) = true;

%% Cargas volumétricas
R = zeros(nNod,nDofNod);
for iele = 1:nel
    nodesEle = nodes(elements(iele,:),:);
    for ipg = 1:npg
        % Punto de Gauss
        ksi = upg(ipg,1);
        eta = upg(ipg,2);  
        % Derivadas de las funciones de forma respecto de ksi, eta
        N = 1/4*[(1 - ksi)*(1 - eta) (1 + ksi)*(1 - eta) (1 + ksi)*(1 + eta) (1 - ksi)*(1 + eta)];
        % Derivadas de las funciones de forma respecto de ksi, eta
        dN = 1/4*[-(1-eta)   1-eta    1+eta  -(1+eta)
                  -(1-ksi) -(1+ksi)   1+ksi    1-ksi ];  
        % Derivadas de x,y, respecto de ksi, eta
        jac = dN*nodesEle;                      
        R(elements(iele,:),2) = R(elements(iele,:),2) + N'*rho*g*t*wpg(ipg)*det(jac);
    end
end

%% Cargas superficiales
a   = 1/sqrt(3);
% Ubicaciones puntos de Gauss
upg = [-a  a];    
% Número de puntos de Gauss
npg = size(upg,2);
wpg = ones(npg,1);
q1 = -480e3; % [N/m^2]
q2 =@(y) -550e3*(1-y/52.5); % [N/m^2]
Q1 = q1*15;
Q2 = integral(q2,0,52.5);
qcheck1 = 0;
qcheck2 = 0;
for iele = 1:nel
    nodesEle = nodes(elements(iele,:),:);
    if nodesEle(4,2)==60
        for ipg = 1:npg
            % Punto de Gauss
            ksi = upg(ipg);
            eta = 1;  
            dN = 1/4*[-(1-eta)   1-eta    1+eta  -(1+eta)
                      -(1-ksi) -(1+ksi)   1+ksi    1-ksi ];  
            % Derivadas de x,y, respecto de ksi, eta
            jac = dN*nodesEle;
            N = 1/4*[(1 - ksi)*(1 - eta) (1 + ksi)*(1 - eta) (1 + ksi)*(1 + eta) (1 - ksi)*(1 + eta)];
            R(elements(iele,4),2) = R(elements(iele,4),2) + N(4)*t*wpg(ipg)*q1*jac(1,1);
            R(elements(iele,3),2) = R(elements(iele,3),2) + N(3)*t*wpg(ipg)*q1*jac(1,1);
            qcheck1 = qcheck1 + (N(4)+N(3))*t*wpg(ipg)*q1*jac(1,1);
        end
    end
%     if nodesEle(3,1)==30&nodesEle(3,2)<=52.5
%         for ipg = 1:npg
%             % Punto de Gauss
%             ksi = 1;
%             eta = upg(ipg);  
%             % Derivadas de las funciones de forma respecto de ksi, eta
%             N = 1/4*[(1 - ksi)*(1 - eta) (1 + ksi)*(1 - eta) (1 + ksi)*(1 + eta) (1 - ksi)*(1 + eta)];
%             % Derivadas de las funciones de forma respecto de ksi, eta
%             dN = 1/4*[-(1-eta)   1-eta    1+eta  -(1+eta)
%                       -(1-ksi) -(1+ksi)   1+ksi    1-ksi ];  
%             % Derivadas de x,y, respecto de ksi, eta
%             delta = (nodesEle(3,2)-nodesEle(2,2))/2;
%             centronodo = (nodesEle(3,2)+nodesEle(2,2))/2;
%             jac = dN*nodesEle;                      
%             R(elements(iele,2),1) = R(elements(iele,2),1) + N(2)*t*wpg(ipg)*q2(centronodo+delta*eta)*jac(2,2);
%             R(elements(iele,3),1) = R(elements(iele,3),1) + N(3)*t*wpg(ipg)*q2(centronodo+delta*eta)*jac(2,2);
%             qcheck2 = qcheck2 + (N(2)+N(3))*t*wpg(ipg)*q2(centronodo+delta*eta)*jac(2,2);
%         end
%     end
    if nodesEle(3,1)==30&nodesEle(3,2)<=52.5
        for ipg = 1:npg
            % Punto de Gauss
            ksi = 1;
            eta = upg(ipg);  
            % Derivadas de las funciones de forma respecto de ksi, eta
            N = 1/4*[(1 - ksi)*(1 - eta) (1 + ksi)*(1 - eta) (1 + ksi)*(1 + eta) (1 - ksi)*(1 + eta)];
            % Derivadas de las funciones de forma respecto de ksi, eta
            dN = 1/4*[-(1-eta)   1-eta    1+eta  -(1+eta)
                      -(1-ksi) -(1+ksi)   1+ksi    1-ksi ];  
            % Derivadas de x,y, respecto de ksi, eta
            jac = dN*nodesEle;
            N = [N(2) N(3)];
            Q = [q2(nodesEle(2,2)) q2(nodesEle(3,2))]';
            R(elements(iele,2:3),1) = R(elements(iele,2:3),1) + N'*N*t*wpg(ipg)*Q*jac(2,2);
            qcheck2 = qcheck2 + sum(N'*N*t*wpg(ipg)*Q*jac(2,2));
        end
    end
end
%% Reduccion Matriz
isFixed = reshape(bc',[],1);
isFree = ~isFixed;

Rr = reshape(R',[],1);

%% Solver
Dr = K(isFree,isFree)\Rr(isFree);

% Reconstrución
D = zeros(nDofTot,1);
D(isFree) = D(isFree) + Dr;

meshplot(elements,nodes+reshape(D,nDofNod,nNod)'*1000,'r')
Dr = reshape(D,nDofNod,nNod)';

uMax = max(abs(Dr(:,1)));
nodoUMax = find( abs(Dr(:,1)) == uMax)
uMax = Dr(nodoUMax,1)
vMax = max(abs(Dr(:,2)));
nodoVMax = find( abs(Dr(:,2)) == vMax)
vMax = Dr(nodoVMax,2)