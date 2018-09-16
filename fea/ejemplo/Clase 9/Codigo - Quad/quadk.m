% ELEMENTOS QUADRILÁTEROS
clear
clc
close all
format short g

eleType = 'Q9'; %'Q4', 'Q8', 'Q9'

% Discretizacion
nodes = [-1.0 -1.0
          1.0 -1.0
          1.0  1.0
         -1.0  1.0
          0.0 -1.0
          1.0  0.0
          0.0  1.0
         -1.0  0.0
          0.0  0.0];

% Q4
%elements = [1 2 3 4];   % Matriz de conectividades - ajustar de acuerdo a número de nodos por elemento 
% Q8
% elements = [1 2 3 4 5 6 7 8];   % Matriz de conectividades - ajustar de acuerdo a número de nodos por elemento
% Q9
elements = [1 2 3 4 5 6 7 8 9];   % Matriz de conectividades - ajustar de acuerdo a número de nodos por elemento


nDofNod = 2;                    % Número de grados de libertad por nodo
nel = size(elements,1);         % Número de elementos
nNod = size(nodes,1);           % Número de nodos
nNodEle = size(elements,2);     % Número de nodos por elemento
nDofTot = nDofNod*nNod;         % Número de grados de libertad
nDims = size(nodes,2);          % Número de dimensiones del problema

% Propiedades del Material
t = 1; % Espesor
E = 1;
NU = 0.3;
C = E/(1 - NU^2)*[ 1.0     NU         0.0
                    NU    1.0         0.0
                   0.0    0.0     (1 - NU)/2 ];
  
              
rsInt = ones(1,2)*2; %*1,*2,*3
[wpg, upg, npg] = gauss(rsInt);

% Matriz de rigidez
K = zeros(nDofTot);
for iele = 1:nel
    Ke = zeros(nDofNod*nNodEle);
    nodesEle = nodes(elements(iele,:),:);
    for ipg = 1:npg
        % Punto de Gauss
        ksi = upg(ipg,1);
        eta = upg(ipg,2);  
        % Derivadas de las funciones de forma respecto de ksi, eta
        dN = shapefunsder([ksi eta],eleType);  
        % Derivadas de x,y, respecto de ksi, eta
        jac = dN*nodesEle;                      
        % Derivadas de las funciones de forma respecto de x,y.
        dNxy = jac\dN;          % dNxy = inv(jac)*dN
        
        B = zeros(size(C,2),nDofNod*nNodEle);
        B(1,1:2:2*nNodEle-1) = dNxy(1,:);
        B(2,2:2:2*nNodEle) = dNxy(2,:);
        B(3,1:2:2*nNodEle-1) = dNxy(2,:);
        B(3,2:2:2*nNodEle) = dNxy(1,:); 

        Ke = Ke + t*B'*C*B*wpg(ipg)*det(jac);
    end
    eleDofs = node2dof(elements(iele,:),nDofNod);
    K(eleDofs,eleDofs) = K(eleDofs,eleDofs) + Ke;  
end