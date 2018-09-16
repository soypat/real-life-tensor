clear
close all
format short g
clc

% Discretizacion
nodes = [ 0.0  0.0
          1.0  0.0
          2.0  0.0
          3.0  0.0
          0.0  1.5
          1.25 0.75
          2.0  1.0
          3.0  1.25
          0.0  3.5
          1.0  3.0
          2.5  2.5
          3.0  3.0
          0.0  5.0
          1.5  5.0
          2.0  5.0
          3.0  5.0];

elements = [1  2  6  5
            2  3  7  6
            3  4  8  7
            5  6 10  9
            6  7 11 10
            7  8 12 11
            9 10 14 13
           10 11 15 14
           11 12 16 15];

nDofNod = 2;                    % Número de grados de libertad por nodo
nNodEle = 4;                    % Número de nodos por elemento
nel = size(elements,1);         % Número de elementos
nNod = size(nodes,1);           % Número de nodos
nDofTot = nDofNod*nNod;         % Número de grados de libertad

bc = false(nNod,nDofNod);       % Matriz de condiciones de borde
bc(1,1:2) = true;
bc(2:4,2) = true;

R = zeros(nNod,nDofNod);        % Vector de cargas
R(13,1)=(nodes(14,1)/2-nodes(13,1)/2);
R(14,1)=(nodes(14,1)/2-nodes(13,1)/2)+(nodes(15,1)/2-nodes(14,1)/2);
R(15,1)=(nodes(15,1)/2-nodes(14,1)/2)+(nodes(16,1)/2-nodes(15,1)/2);
R(16,1)=(nodes(16,1)/2-nodes(15,1)/2);
R(1,1)=-(nodes(2,1)/2-nodes(1,1)/2);
R(2,1)=-(nodes(2,1)/2-nodes(1,1)/2)-(nodes(3,1)/2-nodes(2,1)/2);
R(3,1)=-(nodes(3,1)/2-nodes(2,1)/2)-(nodes(4,1)/2-nodes(3,1)/2);
R(4,1)=-(nodes(4,1)/2-nodes(3,1)/2);
R(1,2)=-(nodes(5,2)/2-nodes(1,2)/2);
R(5,2)=-(nodes(5,2)/2-nodes(1,2)/2)-(nodes(9,2)/2-nodes(5,2)/2);
R(9,2)=-(nodes(9,2)/2-nodes(5,2)/2)-(nodes(13,2)/2-nodes(9,2)/2);
R(13,2)=-(nodes(13,2)/2-nodes(9,2)/2)+(nodes(14,1)-nodes(13,1));
R(4,2)=(nodes(8,2)/2-nodes(4,2)/2);
R(8,2)=(nodes(8,2)/2-nodes(4,2)/2)+(nodes(12,2)/2-nodes(8,2)/2);
R(12,2)=(nodes(12,2)/2-nodes(8,2)/2)+(nodes(16,2)/2-nodes(12,2)/2);
R(16,2)=(nodes(16,2)/2-nodes(12,2)/2)+(nodes(16,1)-nodes(15,1));
R(14,2)=(nodes(14,1)-nodes(13,1))+(nodes(15,1)-nodes(14,1));
R(15,2)=(nodes(15,1)-nodes(14,1))+(nodes(16,1)-nodes(15,1));

% Propiedades del material
E = 10;
NU = -0.1;

meshplot(elements,nodes,'b')

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
% Número de puntos de Gauss
npg = size(upg,1);
wpg = ones(npg,1);

%% Matriz de rigidez
K = zeros(nDofTot);
DetJac = zeros(nel,nDofNod);
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
        
        DetJac(iele,ipg) = det(jac);
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

% Reconstrucción
D = zeros(nDofTot,1);
D(isFree) = D(isFree) + Dr;

% Reacciones
% Rv = K(isFixed,isFree)*D(isFree);
% reacciones = nan(nDofTot,1);
% reacciones(isFixed) = Rv;
% reacciones = (reshape(reacciones,nDofNod,[]))';

%% Recuperación de tensiones en los nodos
stress = zeros(nel,nNodEle,3);
uNod = [ -1 -1
          1 -1
          1  1
         -1  1 ];
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

%% Configuracion deformada
DE = D;
D = (reshape(D,nDofNod,[]))';
nodePosition = nodes + D(:,1:2);

%Graficación
figure(1)
bandplot(elements,nodePosition,stress(:,:,2),[],'k');
meshplot(elements,nodes,'b')
figure(2)
meshplot(elements,nodes,'b')

