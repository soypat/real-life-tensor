% ELEMENTOS FINITOS I - SEGUNDO PARCIAL - 2016

% Funciones de forma Q4 (numeración como figura 6.2-1 pág 205)
N1 = 1/4*(1 - ksi)*(1 - eta);
N2 = 1/4*(1 + ksi)*(1 - eta);
N3 = 1/4*(1 + ksi)*(1 + eta);
N4 = 1/4*(1 - ksi)*(1 + eta);


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
%% Gauss para regla de 3x3 (numeración de nodos como figura 6.3-3 pág 212)          
a   = sqrt(0.6);
% Ubicaciones puntos de Gauss
upg = [ -a  -a
        -a   0
        -a   a
         0  -a    
         0   0
         0   a 
         a  -a
         a   0
         a   a ];
% Número de puntos de Gauss
npg = size(upg,1);
wpg = [5/9, 5/9, 5/9, 5/9, 8/9, 5/9, 5/9, 5/9, 5/9];

%% PROGRAMA DE ELEMENTOS FINITOS PARA ISOPARAMÉTRICOS Q4

elements = load('Nombre de archivo');

nodes = load('Nombre de archivo');

nDofNod = 2;                    % Número de grados de libertad por nodo
nNodEle = 4;                    % Número de nodos por elemento
nel = size(elements,1);         % Número de elementos
nNod = size(nodes,1);           % Número de nodos
nDofTot = nDofNod*nNod;         % Número de grados de libertad

bc = false(nNod,nDofNod);       % Matriz de condiciones de borde

R = zeros(nNod,nDofNod);        % Vector de cargas

% Propiedades del material
E = ;
NU = ;

meshplot(elements,nodes,'b')

%% Matriz Constitutiva (plane strain)
 
% C_Strain = E/((1 + NU)*(1 - 2*NU))*[ 1 - NU      NU            0.0
%                                       NU       1 - NU         0.0
%                                       0.0        0.0     (1 - 2*NU)/2 ];
                          
% Matriz Constitutiva (plane stress)

% C_Stress = E/(1 - NU^2)*[ 1.0     NU         0.0
%                           NU      1.0        0.0
%                           0.0     0.0     (1 - NU)/2 ];
              
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

%% Reduccion Matriz
isFixed = reshape(bc',[],1);
isFree = ~isFixed;

Rr = reshape(R',[],1);

% Solver
Dr = K(isFree,isFree)\Rr(isFree);

% Reconstrución
D = zeros(nDofTot,1);
D(isFree) = D(isFree) + Dr;

%% Configuracion deformada
D = (reshape(D,nDofNod,[]))';
nodePosition = nodes + D(:,1:2);

%Gráfico
bandplot(elements,nodePosition,stress(:,:,2),[],'k');
meshplot(elements,nodes,'b')
