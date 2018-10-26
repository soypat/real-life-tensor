%% An application with high stress gradient
close all; clear; clc
%% Datos
% 1.Acero 2.Aluminio
E = [200 70]*1e3; % [MPa]
NU = [.29 .33];
alfa = [12 24]*1e-6; %[1/C]
eleT = 'Q8';
dT = 1000; % [C]
t = 20; % [mm]
%% Nodos y elementos

elements = load('elembi.txt');
elements = elements(:,2:9);

nodes = load('nodosbi.txt'); % [mm]
nodes = nodes(:,2:3);

nDofNod = 2;                    % Número de grados de libertad por nodo
nNodEle = 8;                    % Número de nodos por elemento
nel = size(elements,1);         % Número de elementos
nNod = size(nodes,1);           % Número de nodos
nDofTot = nDofNod*nNod;         % Número de grados de libertad

bc = false(nNod,nDofNod);       % Matriz de condiciones de borde
R = zeros(nNod,nDofNod);        % Vector de cargas

%meshplot(elements,nodes,'b')

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

%% Matriz de rigidez y deformaciones térmicas
K = zeros(nDofTot);
nodeDofs = reshape(1:nDofTot,nDofNod,nNod)';
epsilon0 = dT*alfa;
P = zeros(nDofTot,1);
for iele = 1:nel
    Ke = zeros(nDofNod*nNodEle);
    Pe = zeros(nDofNod*nNodEle,1);
    nodesEle = nodes(elements(iele,:),:);
    if nodesEle(4,2) == 60
        C = E(1)/(1 - NU(1)^2)*[ 1.0     NU(1)       0.0
                                NU(1)     1.0        0.0
                                 0.0      0.0    (1 - NU(1))/2 ];
        eps = epsilon0(1);
        epsilon = [eps;eps;0];
    else
        C = E(2)/(1 - NU(2)^2)*[ 1.0     NU(2)       0.0
                                NU(2)     1.0        0.0
                                 0.0      0.0    (1 - NU(2))/2 ];
        eps = epsilon0(2);
        epsilon = [eps;eps;0];
    end
    for ipg = 1:npg
        ksi = upg(ipg,1);
        eta = upg(ipg,2);  
        dN = shapefunsder([ksi eta],eleT); 
        jac = dN*nodesEle;                      
        dNxy = jac\dN;
        
        B = zeros(size(C,2),nDofNod*nNodEle);
        B(1,1:2:15) = dNxy(1,:);
        B(2,2:2:16) = dNxy(2,:);
        B(3,1:2:15) = dNxy(2,:);
        B(3,2:2:16) = dNxy(1,:); 

        Ke = Ke + B'*C*B*wpg(ipg)*det(jac)*t;
        Pe = Pe + B'*C*epsilon*wpg(ipg)*det(jac)*t;
    end
    eleDofs = nodeDofs(elements(iele,:),:);
    eleDofs = reshape(eleDofs',[],1);
    K(eleDofs,eleDofs) = K(eleDofs,eleDofs) + Ke;
    P(eleDofs) = P(eleDofs) + Pe;
end


%% BC
bc(1,:) = true;
bc(2,2) = true;

%% Reduccion Matriz
isFixed = reshape(bc',[],1);
isFree = ~isFixed;

Rr = reshape(P,[],1);

% Solver
Dr = K(isFree,isFree)\Rr(isFree);

% Reconstrución
D = zeros(nDofTot,1);
D(isFree) = D(isFree) + Dr;

%% Tensiones en los nodos
stress = zeros(nel,nNodEle,3);
unod = [ -1 -1
          1 -1
          1  1
         -1  1
          0 -1
          1  0
          0  1
         -1  0];
for iele = 1:nel
    nodesEle = nodes(elements(iele,:),:);
    if nodesEle(4,2) == 60
        C = E(1)/(1 - NU(1)^2)*[ 1.0     NU(1)       0.0
                                NU(1)     1.0        0.0
                                 0.0      0.0    (1 - NU(1))/2 ];
        eps = epsilon0(1);
        epsilon = [eps;eps;0];
    else
        C = E(2)/(1 - NU(2)^2)*[ 1.0     NU(2)       0.0
                                NU(2)     1.0        0.0
                                 0.0      0.0    (1 - NU(2))/2 ];
        eps = epsilon0(2);
        epsilon = [eps;eps;0];
    end
    for in = 1:nNodEle
        % Punto de Gauss
        ksi = unod(in,1);
        eta = unod(in,2);  
        % Derivadas de las funciones de forma respecto de ksi, eta
        dNke = shapefunsder([ksi eta],eleT);
        % Derivadas de x,y, respecto de ksi, eta
        J = dNke*nodesEle;                      
        % Derivadas de las funciones de forma respecto de x,y.
        dNxy = J\dNke;
        B = zeros(size(C,2),nDofNod*nNodEle);
        B(1,1:2:15) = dNxy(1,:);
        B(2,2:2:16) = dNxy(2,:);
        B(3,1:2:15) = dNxy(2,:);
        B(3,2:2:16) = dNxy(1,:);
        dofs = reshape(nodeDofs(elements(iele,:),:)',[],1);
        stress(iele,in,:) = C*(B*D(dofs)-epsilon);
    end
end

%% Tensiones en los puntos de superconvergencia
RSextrapolation = unod*sqrt(3);
a   = 1/sqrt(3);
% Ubicaciones puntos de Gauss
upg = [ -a  -a
        -a   a
         a  -a
         a   a ];    
% Número de puntos de Gauss
npg = size(upg,1);
wpg = ones(npg,1);
stress2 = zeros(nel,nNodEle,3);
for iele = 1:nel
    nodesEle = nodes(elements(iele,:),:);
    stressgauss = zeros(npg,3);
    if nodesEle(4,2) == 60
        C = E(1)/(1 - NU(1)^2)*[ 1.0     NU(1)       0.0
                                NU(1)     1.0        0.0
                                 0.0      0.0    (1 - NU(1))/2 ];
        eps = epsilon0(1);
        epsilon = [eps;eps;0];
    else
        C = E(2)/(1 - NU(2)^2)*[ 1.0     NU(2)       0.0
                                NU(2)     1.0        0.0
                                 0.0      0.0    (1 - NU(2))/2 ];
        eps = epsilon0(2);
        epsilon = [eps;eps;0];
    end
    for ipg = 1:npg
        % Punto de Gauss
        ksi = upg(ipg,1);
        eta = upg(ipg,2);  
        % Derivadas de las funciones de forma respecto de ksi, eta
        dNke = shapefunsder([ksi eta],eleT);
        % Derivadas de x,y, respecto de ksi, eta
        J = dNke*nodesEle;                      
        % Derivadas de las funciones de forma respecto de x,y.
        dNxy = J\dNke;
        B = zeros(size(C,2),nDofNod*nNodEle);
        B(1,1:2:15) = dNxy(1,:);
        B(2,2:2:16) = dNxy(2,:);
        B(3,1:2:15) = dNxy(2,:);
        B(3,2:2:16) = dNxy(1,:);
        dofs = reshape(nodeDofs(elements(iele,:),:)',[],1);
        stressgauss(ipg,:) = (C*(B*D(dofs)-epsilon))';
    end
    StressElem = zeros(nNodEle,3);
    for iNod = 1:8
        r = RSextrapolation(iNod,1);
        s = RSextrapolation(iNod,2);
        N = shapefuns([r s],'Q4');
        StressElem(iNod,:) = N*stressgauss;
    end
    stress2(iele,:,:) = StressElem;
end

%% Configuracion deformada
D = (reshape(D,nDofNod,[]))';
nodePosition = nodes + D(:,1:2)*10;

%Gráfico
bandplot(elements,nodePosition,stress2(:,:,1),[],'k');
meshplot(elements,nodes,'b')


































