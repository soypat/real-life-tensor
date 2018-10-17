clear
clc
close all

eleType = 'Q4';% 'Q4', 'Q9'
escala=1; %Escala desplazamientos
%% Discretizacion

nodos=load('NodosT1.txt');


elementos=load('ElementosT1.txt');
% elementos = 1:4;   % Matriz de conectividades - ajustar de acuerdo a número de nodos por elemento
% elementos=[1 5 2 6 3 7 4 8];
%% Propiedades del Material
E=1;
nu=0.3;
t=1; %espesor
%Plane Strain
C = (E/((1+nu)*(1-2*nu)))*[1-nu    nu      0;
                           nu  1-nu      0;
                            0    0  0.5-nu];
%Plane Stress
% C = (E/(1-nu^2))*[1 nu 0;
%                   nu 1 0;
%                   0 0 (1-nu)/2];


%% Definiciones
Ndofpornod = 2;                    % grados de libertad por nodo
Nelem = size(elementos,1);         % elementos
Nnod = size(nodos,1);           % nodos
Nnodporelem = size(elementos,2);     % nodos por elemento
doftot = Ndofpornod*Nnod;         % grados de libertad
Ndims = 2;          % dimensiones del problema

%% Puntos de Gauss
rsInt = 2*ones(1,2);
[wpg, upg, npg] = gauss(rsInt);

%% Matriz de rigidez
K = zeros(doftot);
for iele = 1:Nelem
    Ke = zeros(Ndofpornod*Nnodporelem);
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
        
        B = zeros(size(C,2),Ndofpornod*Nnodporelem);
        B(1,1:2:Ndofpornod*Nnodporelem-1) = dNxy(1,:);
        B(2,2:2:Ndofpornod*Nnodporelem) = dNxy(2,:);
        B(3,1:2:Ndofpornod*Nnodporelem-1) = dNxy(2,:);
        B(3,2:2:Ndofpornod*Nnodporelem) = dNxy(1,:);
        
        Djac = det(jac);
        Ke = Ke + B'*C*B*wpg(ipg)*Djac*t;
    end
    eleDofs = node2dof(elementos(iele,:),Ndofpornod);
    K(eleDofs,eleDofs) = K(eleDofs,eleDofs) + Ke;
end

%% Condiciones de borde
bc = false(Nnod,Ndofpornod);       % Matriz de condiciones de borde
% bc(1,1:2) = true;
% bc(2,2) = true;
for inod=1:Nnod
    if nodos(inod,2)==0
        bc(inod,[1 2])=true;
    end
end

%%
figure(1)%Meshplot(elementos,nodos,bc,color,verbatim)
Meshplot(elementos,nodos,bc,'k',1)
% MeshplotTrigMec(elementos,nodos,bc,'k',1)

%% Cargas de linea
R = zeros(doftot,1);% Vector de cargas
DOF=reshape(1:doftot,2,[])';
surfnod=[4 3];
wpg=[1 1];
a=-sqrt(3)^-1;
upg=[-a 1;a 1];
npg=2;
sig=-.02;
for e=1:Nelem
    index=elementos(e,:);
    meindof=DOF(index,:);
    elenod=nodos(index,:);
    A=elenod(:,2)==60; % y=60 => Cargo
    if  isempty(find(A,1))
        continue
    end
    for ipg=1:npg
        ksi = upg(ipg,1);
        eta = upg(ipg,2);
        N = shapefuns([ksi eta],eleType);
        dN = shapefunsder([ksi eta],eleType);
        jac = dN*elenod;
        W=wpg(ipg);
        for snod=surfnod
            R(meindof(snod,1))=R(meindof(snod,1))+N(snod)*(-jac(1,2))*sig*W;
            R(meindof(snod,2))=R(meindof(snod,2))+N(snod)*jac(1,1)*sig*W;
        end
    end
end
R=reshape(R,2,[])';
% Reduccion Matriz
isFixed = reshape(bc',[],1);
isFree = ~isFixed;

Rr = reshape(R',[],1);

% Solver
Dr = K(isFree,isFree)\Rr(isFree);

% Reconstrucción
D = zeros(doftot,1);
D(isFree) = D(isFree) + Dr;

% Reacciones
Rv = K(isFixed,isFree)*D(isFree);
reacciones = nan(doftot,1);
reacciones(isFixed) = Rv;
reacciones = (reshape(reacciones,Ndofpornod,[]))';

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
stress = zeros(Nnodporelem,3,Nelem);
for iele = 1:Nelem
    nodesEle = nodos(elementos(iele,:),:);
    index=elementos(iele,:);
    meindof=DOF(index,:);
    elenod=nodos(index,:);
    eleDofs = reshape(meindof',[],1);
%     eleDofs = node2dof(elementos(iele,:),Ndofpornod);
    u=1;
    for inode = 1:Nnodporelem
        % Punto de Gauss
        ksi = uNod(inode,1);
        eta = uNod(inode,2);
        % Derivadas de las funciones de forma respecto de ksi, eta
        dN  = shapefunsder([ksi eta],eleType);
        % Derivadas de x,y, respecto d % '2'; % e ksi, eta
        jac  = dN*nodesEle;
        % Derivadas de las funciones de forma respecto de x,y.
        dNxy = jac\dN;          % dNxy = inv(jac)*dN
        
        B = zeros(size(C,2),Ndofpornod*Nnodporelem);
        B(1,1:2:Ndofpornod*Nnodporelem-1) = dNxy(1,:);
        B(2,2:2:Ndofpornod*Nnodporelem) = dNxy(2,:);
        B(3,1:2:Ndofpornod*Nnodporelem-1) = dNxy(2,:);
        B(3,2:2:Ndofpornod*Nnodporelem) = dNxy(1,:);
        
        
        stress(inode,:,iele) = C*B*D(eleDofs);
    end
end
    
%% Configuración deformada
nodePosition = nodos + escala*(reshape(D,Ndofpornod,[]))';
figure(1)
Meshplot(elementos,nodePosition,bc,'r',0)
%% Gráfico
figure(2)
S=2;
plotvar=squeeze(stress(:,S,:))';
bandplot(elementos,nodePosition,plotvar,[],'k');