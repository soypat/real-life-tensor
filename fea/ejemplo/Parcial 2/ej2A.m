clc
clear

%% Discretizacion 

nodes = load('NodosT1.txt')*1e3;                        % Paso todo a mm.
elements = load('ElementosT1.txt');

meshplot(elements,nodes,'b')

% Trabajo en mm, MPa y N.

nDofNod = 2;                    % Número de grados de libertad por nodo
nNodEle = 4;                    % Número de nodos por elemento
nel = size(elements,1);         % Número de elementos
nNod = size(nodes,1);           % Número de nodos
nDofTot = nDofNod*nNod;         % Número de grados de libertad

bc = false(nNod,nDofNod);       % Matriz de condiciones de borde

R = zeros(nNod,nDofNod);        % Vector de cargas

% Propiedades del material
E = 30e3;
NU = 0.18;
t = 1; %% Ancho unitario
gama = -2e-05;


%% Matriz Constitutiva (plane strain) represa
 
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
    re = zeros(4,1);
    eleNods = elements(iele,:);
    for ipg = 1:npg
        % Punto de Gauss
        ksi = upg(ipg,1);
        eta = upg(ipg,2);  
        % Derivadas de las funciones de forma respecto de ksi, eta
        
        N1 = 1/4*(1 - ksi)*(1 - eta);
        N2 = 1/4*(1 + ksi)*(1 - eta);
        N3 = 1/4*(1 + ksi)*(1 + eta);
        N4 = 1/4*(1 - ksi)*(1 + eta);

        N = [N1,N2,N3,N4];
        
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
        
        
        %  Cargo las Volumetricas aprovechando el ciclo
        
        re = re + N'*N*det(jac)*ones(4,1)*gama*t;
        R(eleNods,2) = re + R(eleNods,2);
        
                
    end
    eleDofs = nodeDofs(elements(iele,:),:);
    eleDofs = reshape(eleDofs',[],1);
    K(eleDofs,eleDofs) = K(eleDofs,eleDofs) + Ke;  
end

movRigidos = find(eig(K)<1e-9);

%% Condiciones de borde y cargas

% Encontramos los nodos que estan fijos

for iNod = 1:nNod
    if abs(nodes(iNod,2)) < 1e-9
        bc(iNod,:) = true;
    end
end

%% Cargas

gy = -0.480;
q = @(y) +(0.55/52.5e3)*y - 0.55;

for iEle = 1:nel
    eleNods = [elements(iEle,:),elements(iEle,1)];
    
    % Superficial
    
    for iNod = 1:4
        if abs(nodes(eleNods(iNod),1)-30e3) <1e-9 && nodes(eleNods(iNod),2)<53e3
            if abs(nodes(eleNods(iNod+1),1)-30e3) <1e-9 && nodes(eleNods(iNod+1),2)<53e3
                nodA = eleNods(iNod);
                nodB = eleNods(iNod + 1);
                
                if nodes(nodA,2)<nodes(nodB,2)
                    a = nodes(nodA,2);
                    b = nodes(nodB,2);
                else
                    a = nodes(nodB,2);
                    b = nodes(nodA,2);
                    k = nodB;
                    nodB = nodA;
                    nodA = k;
                end
                % Simpson, las funciones de forma son lineales, y la carga
                % lineal, orden 2. Alcanza Simpson.
                l = b-a; 
                R(nodA,1) = l/6 * (q(a) + 2*q((a+b)/2))*t + R(nodA,2);
                R(nodB,1) = l/6 * (q(b) + 2*q((a+b)/2))*t + R(nodB,2);
            end
        end
        
        % Para el lado de arriba
        if abs(nodes(eleNods(iNod),2)-60e3) <1e-9 
            if abs(nodes(eleNods(iNod+1),2)-60e3) <1e-9 
                nodA = eleNods(iNod);
                nodB = eleNods(iNod + 1);
                
                if nodes(nodA,1)<nodes(nodB,1)
                    a = nodes(nodA,1);
                    b = nodes(nodB,1);
                else
                    a = nodes(nodB,1);
                    b = nodes(nodA,1);
                    k = nodB;
                    nodB = nodA;
                    nodA = k;
                end
                % Simpson, las funciones de forma son lineales, y la carga
                % cte, orden 1. Alcanza Simpson.
                l = b-a; 
                R(nodA,2) = l/6 * (gy + 2*gy)*t + R(nodA,2);
                R(nodB,2) = l/6 * (gy + 2*gy)*t + R(nodB,2);
            end
        end
        
    end
end




%% Reduccion Matriz
isFixed = reshape(bc',[],1);
isFree = ~isFixed;

Rr = reshape(R',[],1);

% Solver
Dr = K(isFree,isFree)\Rr(isFree);

% Reconstrución
D = zeros(nDofTot,1);
D(isFree) = D(isFree) + Dr;

D = (reshape(D,nDofNod,[]))';

%% Desp Maximos b)

uMax = max(abs(D(:,1)))
nodoUMax = find( abs(D(:,1)) == uMax)
vMax = max(abs(D(:,2)));
nodoVMax = find( abs(D(:,2)) == vMax)
vMax = D(nodoVMax,2)

nodePosition = nodes + D(:,1:2)*100;
meshplot(elements,nodePosition,'r')
        
        
        
                
                
                   
            
        
        
       
      



