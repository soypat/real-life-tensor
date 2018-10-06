%Q8, tau_xy
TestType ='Tau_xy';% 'Sigma_y' ; % 'Sigma_y' ; % 'Sigma_x' ; %
eleType = 'Q8';% 'Q4', 'Q9'
escala=1; %Escala desplazamientos
%% Discretizacion
nodos = [ 0      0.00
          1      0.00
          1       1
          0       1
          .5      0
          1      .5
          .5      1.00
          0      .5];
elementos = 1:8;

%% Propiedades del Material {plain stress}
E=1;
nu=0.3;
C = (E/(1-nu^2))*[1 nu 0;
                  nu 1 0;
                  0 0 (1-nu)/2];


%% Definiciones
nDofNod = 2;                    % grados de libertad por nodo
nel = size(elementos,1);         % elementos
nNod = size(nodos,1);           % nodos
nNodEle = size(elementos,2);     % nodos por elemento
nDofTot = nDofNod*nNod;         % grados de libertad
nDims = size(nodos,2);          % dimensiones del problema

%% Condiciones de borde
bc = false(nNod,nDofNod);       % Matriz de condiciones de borde
bc(1,1:2) = true;
bc(2,2) = true;
%%
figure(1)
Meshplot(elementos,nodos,bc,'k',1)

%% Cargas
R = zeros(nNod,nDofNod);        % Vector de cargas

%% CASO TAU_XY
%busco tauxy=2 unitario
S=3;
R(1,1)=1/6;%Q8 loading
R(1,2)=1/6;
R(5,1)=2/3;
R(2,1)=1/6; 
R(2,2)=-1/6;
R(6,2)=-2/3;
R(3,1)=-1/6;
R(3,2)=-1/6;
R(7,1)=-2/3;
R(4,1)=-1/6;
R(4,2)=-1/6;
R(8,2)=-2/3;
R=R*2;
%% Puntos de Gauss
rsInt = 3*ones(1,2); %2 puntos de gauss. Recibe cuantos 
[wpg, upg, npg] = gauss(rsInt);

%% Matriz de rigidez
K = zeros(nDofTot);
A = 0; %
jmin = 1E10;%guardamos j min for our information
for iele = 1:nel
    Ke = zeros(nDofNod*nNodEle);
    nodesEle = nodos(elementos(iele,:),:);
    for ipg = 1:npg %ultimo punto de gauss
        % Punto de Gauss
        ksi = upg(ipg,1);
        eta = upg(ipg,2);
        
        % Funciones de forma respecto de ksi, eta
        N = shapefuns([ksi eta],eleType);
        %shapefuns([-1 1],eleType); nos da [0 0 0 1] porque estamos parado
        %en el lado superior izquierdo del elemento
        % Derivadas de las funciones de forma respecto de ksi, eta
        dN = shapefunsder([ksi eta],eleType);
        %sum(dN,2) siempre da uno
        % Derivadas de x,y, respecto de ksi, eta
        jac = dN*nodesEle; %Mucho lio para UN punto del elemento
        %dN es el mismo para el elemento, lo que cambia es nodesEle para
        %obtener el jabociano
        % Derivadas de las funciones de forma respecto de x,y.
        dNxy = jac\dN;          % dNxy = inv(jac)*dN
        
        B = zeros(size(C,2),nDofNod*nNodEle);
        B(1,1:2:nDofNod*nNodEle-1) = dNxy(1,:);
        B(2,2:2:nDofNod*nNodEle) = dNxy(2,:);
        B(3,1:2:nDofNod*nNodEle-1) = dNxy(2,:);
        B(3,2:2:nDofNod*nNodEle) = dNxy(1,:);
        
        Djac = det(jac);
        Ke = Ke + B'*C*B*wpg(ipg)*Djac;
        A = A + wpg(ipg)*Djac;
        if Djac < jmin
            jmin = Djac;
        end
    end
    eleDofs = node2dof(elementos(iele,:),nDofNod);
    K(eleDofs,eleDofs) = K(eleDofs,eleDofs) + Ke;
end

% Reduccion Matriz
isFixed = reshape(bc',[],1);
isFree = ~isFixed;

Rr = reshape(R',[],1);

% Solver
Dr = K(isFree,isFree)\Rr(isFree);

% Reconstrucción
D = zeros(nDofTot,1);
D(isFree) = D(isFree) + Dr;

% Reacciones
Rv = K(isFixed,isFree)*D(isFree);
reacciones = nan(nDofTot,1);
reacciones(isFixed) = Rv;
reacciones = (reshape(reacciones,nDofNod,[]))';

%% Recuperación de tensiones en los nodos
uNod = [-1 -1
         1 -1
         1  1
        -1  1
         0 -1
         1  0
         0  1
        -1  0];%Q8 uNod

stress = zeros(nNodEle,3,nel);
for iele = 1:nel
    nodesEle = nodos(elementos(iele,:),:);
    for inode = 1:nNodEle
        % Punto de Gauss
        ksi = uNod(inode,1);
        eta = uNod(inode,2);
        % Derivadas de las funciones de forma respecto de ksi, eta
        dN  = shapefunsder([ksi eta],eleType);
        % Derivadas de x,y, respecto d % '2'; % e ksi, eta
        jac  = dN*nodesEle;
        % Derivadas de las funciones de forma respecto de x,y.
        dNxy = jac\dN;          % dNxy = inv(jac)*dN
        
        B = zeros(size(C,2),nDofNod*nNodEle);
        B(1,1:2:nDofNod*nNodEle-1) = dNxy(1,:);
        B(2,2:2:nDofNod*nNodEle) = dNxy(2,:);
        B(3,1:2:nDofNod*nNodEle-1) = dNxy(2,:);
        B(3,2:2:nDofNod*nNodEle) = dNxy(1,:);
        
        eleDofs = node2dof(elementos(iele,:),nDofNod);
        stress(inode,:,iele) = C*B*D(eleDofs);% Aca esta resuelto el problem
        %C*B es sigma, por los desplazaments
    end
end
    
%% Configuración deformada
nodePosition = nodos + escala*(reshape(D,nDofNod,[]))';
figure(1)
Meshplot(elementos,nodePosition,bc,'r',0)
%% Gráfico
figure(2)
bandplot(elementos,nodePosition,stress(:,S,:)',[],'k');