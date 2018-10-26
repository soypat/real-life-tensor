%% Problema 1 P15
clear; close all; clc
format short g
%% Datos
E = 1000;
NU = 0.3;
P =700;

%% Nodos y elementos
load('Nodos5x7Q4');
load('Elementos5x7Q4');
elements = elements(1:end-2,:);
elements = elements([1:26 29:33],:);
nodes = nodes(~(nodes(:,1)>5&nodes(:,2)>3),:);

meshplot(elements,nodes, 'b');

nDofNod = 2;                    % Número de grados de libertad por nodo
nNodEle = 4;                    % Número de nodos por elemento
nel = size(elements,1);         % Número de elementos
nNod = size(nodes,1);           % Número de nodos
nDofTot = nDofNod*nNod;         % Número de grados de libertad

bc = false(nNod,nDofNod);       % Matriz de condiciones de borde
bc(nodes(:,1)-0<=1e-4,1:2) = true;
bc(nodes(:,2)-0<=1e-4,1:2) = true;
bc(7:8,:) = false;

R = zeros(nNod,nDofNod);
R(8,2) = -P;


%% Matriz Constitutiva (plane stress)
C = E/(1 - NU^2)*[ 1.0     NU         0.0
                    NU    1.0         0.0
                   0.0    0.0     (1 - NU)/2 ];

%% Matriz de rigidez
% Gauss
[wpg, upg, npg] = gauss([2 2]);
upg = upg([1 3 4 2],:);
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
        dN = shapefunsder([ksi eta],'Q4');  
        % Derivadas de x,y, respecto de ksi, eta
        jac = dN*nodesEle;                      
        % Derivadas de las funciones de forma respecto de x,y.
        dNxy = jac\dN;          % dNxy = inv(jac)*dN
        
        B = zeros(size(C,2),nDofNod*nNodEle);
        B(1,1:2:nNodEle*2-1) = dNxy(1,:);
        B(2,2:2:nNodEle*2)   = dNxy(2,:);
        B(3,1:2:nNodEle*2-1) = dNxy(2,:);
        B(3,2:2:nNodEle*2)   = dNxy(1,:); 

        Ke = Ke + B'*C*B*wpg(ipg)*det(jac);
    end
    eleDofs = nodeDofs(elements(iele,:),:);
    eleDofs = reshape(eleDofs',[],1);
    K(eleDofs,eleDofs) = K(eleDofs,eleDofs) + Ke;  
end
% Constraint


%% Reduccion Matriz
isFixed = logical(reshape(bc',[],1));
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

%% Graficos
figure
meshplot(elements,nodes,'b')
hold on
Dreshaped = (reshape(D,nDofNod,[]))';
nodePosition = nodes + Dreshaped(:,1:2);
meshplot(elements,nodePosition,'r')