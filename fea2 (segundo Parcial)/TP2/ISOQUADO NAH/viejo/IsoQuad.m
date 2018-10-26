clear
% clc
close all

TestType = 'Sigma_y' ; % 'Sigma_y' ; % 'Sigma_x' ; %
eleType = 'Q8';% 'Q4', 'Q9'
escala=1; %Escala desplazamientos
%% Discretizacion

% CARGAR MALLA de Adina
% if strcmp(eleType,'Q4')
%     nodos = [ 0.00      0.00
%               1.00      0.00
%               1.00      1.00
%               0.00      1.00];
%     elementos = 1:4;   % Matriz de conectividades - ajustar de acuerdo a número de nodos por elemento
% elseif strcmp(eleType,'Q8')
%     nodos = [ 0.00      0.00
%               1.00      0.00
%               1.00      1.00
%               0.00      1.00
%               0.50      0.00
%               1.00      0.50
%               0.50      1.00
%               0.00      0.50];
%     elementos = 1:8;
% elseif strcmp(eleType,'Q9')
%     nodos = [ 0.00      0.00
%               1.00      0.00
%               1.00      1.00
%               0.00      1.00
%               0.50      0.00
%               1.00      0.50
%               0.50      1.00
%               0.00      0.50
%               0.50      0.50];
%     elementos = 1:9;
% end

%% Propiedades del Material
E=1;
nu=0.3;
%Plane Strain
%C = (E/((1+nu)*(1-2*nu)))*[1-nu    nu      0;
                           %nu  1-nu      0;
                            %0    0  0.5-nu];
%Plane Stress
C = (E/(1-nu^2))*[1 nu 0;
                  nu 1 0;
                  0 0 (1-nu)/2];


%% Definiciones
%reescribimos la matriz elementos y nodos, sino lo cargo de adina
elementos=[1 5 9 8;
    5 2 6 9;
    9 6 3 7;
    8 9 7 4];

nodos = [ 0.00      0.00
    1.00      0.00
    1.00      1.00
    0.00      1.00
    0.50      0.00
    1.00      0.50
    0.50      1.00
    0.00      0.50
    0.50      0.50];
       
       
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
% R = zeros(nNod,nDofNod);        % Vector de cargas
% switch TestType
%     case 'Sigma_x' %sigma x constante
%         S=1;
%         if strcmp(eleType,'Q4')
%             R(1,1)=-1/2;
%             R(2,1)=1/2;
%             R(3,1)=1/2;
%             R(4,1)=-1/2;
%         else
%             R(1,1) = -1/6;
%             R(8,1) = -2/3;
%             R(4,1) = -1/6;
%             R(2,1) = 1/6;
%             R(6,1) = 2/3;
%             R(3,1) = 1/6;
%         end
%         
%     case 'Sigma_y' %sigma y constante
%         S=2;
%         if strcmp(eleType,'Q4')
%             R(1,2)=-1/2;
%             R(2,2)=-1/2;
%             R(3,2)=1/2;
%             R(4,2)=1/2;
%         else
%             R(1,2) = -1/6;
%             R(5,2) = -2/3;
%             R(2,2) = -1/6;
%             R(3,2) = 1/6;
%             R(7,2) = 2/3;
%             R(4,2) = 1/6;
%         end
%         
%     case 'Tau_xy' %corte xy constante
%         S=3;
%         if strcmp(eleType,'Q4')
%             R(1,1)=1/2;
%             R(1,2)=1/2;
%             R(2,1)=1/2;
%             R(2,2)=-1/2;
%             R(3,1)=-1/2;
%             R(3,2)=-1/2;
%             R(4,1)=-1/2;
%             R(4,2)=1/2;
%         else
%             R(1,1)=1/6;
%             R(1,2)=1/6;
%             R(5,1)=2/3;
%             R(2,1)=1/6;
%             R(2,2)=-1/6;
%             R(6,2)=-2/3;
%             R(3,1)=-1/6;
%             R(3,2)=-1/6;
%             R(7,1)=-2/3;
%             R(4,1)=-1/6;
%             R(4,2)=-1/6;
%             R(8,2)=-2/3;
%         end
% end
%las cargo a mano!
%SX
R=zeros(9,2);
S=1;
R(1,1) = -1/8;
R(8,1) = -1/4;
R(4,1) = -1/8;
R(2,1) = 1/8;
R(6,1) = 1/4;
R(3,1) = 1/8;

% %SY
% R(1,2) = 1/4;
% R(5,2) = 1/2;
% R(2,2) = 1/4;
% R(3,2) = -1/4;
% R(7,2) = -1/2;
% R(4,2) = -1/4;
% 
% %Txy
% R(1,1)=1/2+R(1,1);
% R(1,2)=1/2+R(1,2);
% R(5,1)=2/3+R(5,1);
% R(2,1)=1/2+R(2,1);
% R(2,2)=-1/2+R(2,2);
% R(6,2)=-1+R(6,2);
% R(3,1)=-1/2+R(3,1);
% R(3,2)=-1/2+R(3,2);
% R(7,1)=-1+R(7,1);
% R(4,1)=-1/2+R(4,1);
% R(4,2)=1/2+R(4,2);
% R(8,2)=1+R(8,2);

%% Puntos de Gauss
rsInt = 3*ones(1,2);
[wpg, upg, npg] = gauss(rsInt);

%% Matriz de rigidez
K = zeros(nDofTot);
A = 0;
jmin = 1E10;
for iele = 1:nel
    Ke = zeros(nDofNod*nNodEle);
    nodesEle = nodos(elementos(iele,:),:);
    for ipg = 1:npg
        % Punto de Gauss
        ksi = upg(ipg,1);
        eta = upg(ipg,2);
        
        % Funciones de forma respecto de ksi, eta
        N = shapefuns([ksi eta],eleType);
        
        % Derivadas de las funciones de forma respecto de ksi, eta
        dN = shapefunsder([ksi eta],eleType);
        % Derivadas de x,y, respecto de ksi, eta
        jac = dN*nodesEle;
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
if strcmp(eleType,'Q4')
    uNod = [-1 -1
             1 -1
             1  1
            -1  1];
elseif strcmp(eleType,'Q8')
    uNod = [-1 -1
             1 -1
             1  1
            -1  1
             0 -1
             1  0
             0  1
            -1  0];
elseif strcmp(eleType,'Q9')    
    uNod = [-1 -1
             1 -1
             1  1
            -1  1
             0 -1
             1  0
             0  1
            -1  0
             0  0];
end
stress = zeros(nNodEle,nel,3);
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
        stress(inode,iele,:) = C*B*D(eleDofs);
    end
end
S=1;
    
%% Configuración deformada
nodePosition = nodos + escala*(reshape(D,nDofNod,[]))';
figure(1)
Meshplot(elementos,nodePosition,bc,'r',0)
%% Gráfico
figure(2)
bandplot(elementos,nodePosition,stress(:,:,S)',[],'k');

