
clear
clc
close all
format short g

nodes = load('error/nodosO.txt');
nodes(:,[1 2 5]) = [];

elements =load('error/elesO.txt');
elements(:,[1 6:10]) = [];

nDofNod = 2;                    % Numero de grados de libertad por nodo
nNodEle = 4;                    % Numero de nodos por elemento
nel = size(elements,1);         % Numero de elementos
nNod = size(nodes,1);           % Numero de nodos
nDofTot = nDofNod*nNod;         % Numero de grados de libertad

bc = false(nNod,nDofNod);       % Matriz de condiciones de borde
bc([13:16 59:61],2) = true;
bc([43:3:49 52:3:58 61],1) = true;

R = zeros(nNod,nDofNod);        % Vector de cargas
R(29,2) = 1/2;
R([30 31 32 41 42],2) = 1;
R(43,2) = 1/2;

% Propiedades del material
E = 1;
NU = 0.3;

% meshplot(elements,nodes,'b')

%% Matriz Constitutiva (plane stress)

C = E/(1 - NU^2)*[ 1.0     NU         0.0
                    NU    1.0         0.0
                   0.0    0.0     (1 - NU)/2 ];
               

%% Gauss           
a   = 1/sqrt(3);
% Ubicaciones puntos de Gauss
upg = [ -a  -a
         a  -a
         a   a
        -a   a ];    
% Numero de puntos de Gauss
npg = size(upg,1);
wpg = ones(npg,1);

%% Matriz de rigidez
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

%% Reduccion Matriz
isFixed = reshape(bc',[],1);
isFree = ~isFixed;

Rr = reshape(R',[],1);

% Solver
Dr = K(isFree,isFree)\Rr(isFree);

% Reconstruccion
D = zeros(nDofTot,1);
D(isFree) = D(isFree) + Dr;

% Reacciones
Rv = K(isFixed,isFree)*D(isFree);
reacciones = nan(nDofTot,1);
reacciones(isFixed) = Rv;
reacciones = (reshape(reacciones,nDofNod,[]))';

%% Recuperacion de tensiones en los puntos de Gauss
stress = zeros(nel,nNodEle,3);
uNod = [ -a  -a
          a  -a
          a   a
         -a   a ];
for iele = 1:nel
    nodesEle = nodes(elements(iele,:),:);
    for inode = 1:nNodEle
        % Punto de Gauss
        ksi = uNod(inode,1);
        eta = uNod(inode,2);  
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
        
        eleDofs = nodeDofs(elements(iele,:),:);
        eleDofs = reshape(eleDofs',[],1);
        stress(iele,inode,:) = C*B*D(eleDofs);
    end
end

%% tensiones en el centro de cada elemento
stressCentro = zeros(nel,nNodEle,3);
for iele = 1:nel
    nodesEle = nodes(elements(iele,:),:);
    for inode = 1:nNodEle
        % Punto de Gauss
        ksi = 0;
        eta = 0;  
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
        
        eleDofs = nodeDofs(elements(iele,:),:);
        eleDofs = reshape(eleDofs',[],1);
        stressCentro(iele,inode,:) = C*B*D(eleDofs);
    end
end

%% promediado de tensiones en los nodos
avgStress = zeros(nNod,3);
for inode = 1:nNod
    [I,J] = find(elements == inode);
    nShare = length(I);
    for ishare = 1:nShare
        avgStress(inode,:) = avgStress(inode,:) + squeeze(stressCentro(I(ishare),J(ishare),:))';
    end
    avgStress(inode,:) = avgStress(inode,:) / nShare;
end

%% calculo de eta_el, e2, U2
invC = C\eye(3);
eta_el = zeros(nel,1);
e2_el = zeros(nel,1);
U2_el = zeros(nel,1);
for iele = 1:nel
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
        
        % funciones de forma
        N4 = 0.25*(1 - ksi)*(1 + eta);
        N3 = 0.25*(1 + ksi)*(1 + eta);
        N2 = 0.25*(1 + ksi)*(1 - eta);
        N1 = 0.25*(1 - ksi)*(1 - eta);
        N = [ N1 N2 N3 N4 ];
        eleStress =  squeeze(stress(iele,ipg,:));            % tensiones "directas"
        starStress = ( N * avgStress(elements(iele,:),:) )'; % tensiones mejoradas
        
        e2_el(iele) = e2_el(iele) + (starStress - eleStress)' * ... 
                    invC * (starStress - eleStress) * wpg(ipg) * det(jac);
                
        U2_el(iele) = U2_el(iele) + eleStress' * invC * eleStress * ...
                      wpg(ipg) * det(jac);
    end
    
    eta_el(iele) = sqrt( e2_el(iele) / (e2_el(iele) + U2_el(iele)) );
    
end

etaG = sqrt( sum(e2_el) / (sum(e2_el) + sum(U2_el)) );


%% Configuracion deformada
D = (reshape(D,nDofNod,[]))';
nodePosition = nodes + D(:,1:2);

%Graficacion
limites = [2.5E-1,3.5E-1];

figure(1)
bandplot(elements,nodePosition,stressCentro(:,:,2),limites,'k');
title('Tensiones en el centro de cada elemento.')

figure(2)
scalarbandplot(elements,nodePosition,avgStress(:,2),limites,'k',[],'interp');
title('Tensiones promediadas.')

figure(3)
scalarbandplot(elements,nodePosition,eta_el,[],'k',[],'flat');


