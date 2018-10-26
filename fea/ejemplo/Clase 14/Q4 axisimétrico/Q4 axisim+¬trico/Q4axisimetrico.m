%% Modelo Q4 axisimétrico. Arandela girando con velocidad angular constante
clear 
clc
close all

%% Discretización

wt = 10;
Ro = 20;
R1 = Ro + wt;
H = 1;

Nodos = [ Ro 0
          R1 0
          R1 H
          Ro H ];

% Elementos = [1  2  3  4];      %Matriz de conectividades

% load('Elementos10x2.mat','Elementos')
% load('Nodos10x2.mat','Nodos')
% load('Elementos30x6.mat','Elementos')
% load('Nodos30x6.mat','Nodos')
load('Elementos60x10.mat','Elementos')
load('Nodos60x10.mat','Nodos')

eleType = 'Q4';
nDofNod = 2;                    % Número de grados de libertad por nodo
switch eleType
    case 'Q4'
        nNodEle = 4;            % Número de nodos por elemento
end
nel = size(Elementos,1);         % Número de elementos
nNod = size(Nodos,1);           % Número de nodos
nDofTot = nDofNod*nNod;         % Número de grados de libertad

bc = false(nNod,nDofNod);       % Matriz de condiciones de borde
bc(1,2) = true;
% bc(2,2) = true;

w = 1;                          % Velocidad angular
RHO = 1;                        % Densidad

%% Propiedades del material
E = 1000;
NU = 0.33;
%%
% meshplot(Elementos,Nodos,'b')

%% Matriz Constitutiva

f = NU/(1 - NU);
g = (1 - 2*NU)/(2*(1 - NU));
C = [ 1 f f 0
      f 1 f 0
      f f 1 0
      0 0 0 g ] * (1 - NU)*E/((1 + NU)*(1 - 2*NU));

%% Matriz de rigidez
% Gauss           
rsInt = 2*ones(1,2);
[wpg, upg, npg] = gauss(rsInt);
% Funciones de forma y derivadas
Ni = shapefuns   (upg,eleType);
dN = shapefunsder(upg,eleType);

K = zeros(nDofTot);
for iele = 1:nel
    Ke = zeros(nDofNod*nNodEle);
    nodesEle = Nodos(Elementos(iele,:),:);
    for ipg = 1:npg
        % Derivadas de x,y, respecto de ksi, eta
        jac = dN(:,:,ipg)*nodesEle;                      
        % Derivadas de las funciones de forma respecto de x,y.
        dNxy = jac\dN(:,:,ipg);          % dNxy = inv(jac)*dN(:,:,ipg)

        r = Ni(:,:,ipg)*nodesEle(:,1);
        
        B = zeros(size(C,2),nDofNod*nNodEle);
        B(1,1:2:7) = dNxy(1,:);
        B(2,1:2:7) = Ni(:,:,ipg)/r;
        B(3,2:2:8) = dNxy(2,:);
        B(4,1:2:7) = dNxy(2,:);
        B(4,2:2:8) = dNxy(1,:);        
        
        Ke = Ke + B'*C*B*r*det(jac)*wpg(ipg);
    end
    eleDofs = node2dof(Elementos(iele,:),nDofNod);
    K(eleDofs,eleDofs) = K(eleDofs,eleDofs) + Ke;  
end

%% Vector de cargas velocidad angular constante
% Gauss           
rsInt = 2*ones(1,2);
[wpg, upg, npg] = gauss(rsInt);
% Funciones de forma
[Ni N] = shapefuns   (upg,eleType);
   dN  = shapefunsder(upg,eleType);

R = zeros(nDofTot,1);
for iele = 1:nel
    Re = zeros(nDofNod*nNodEle,1);
    nodesEle = Nodos(Elementos(iele,:),:);
    for ipg = 1:npg
        % Derivadas de x,y, respecto de ksi, eta
        jac = dN(:,:,ipg)*nodesEle;                      
        % Radio
        r = Ni(:,:,ipg)*nodesEle(:,1);  
        Re = Re + N(:,:,ipg)'*r^2*RHO*w^2*[1; 0]*det(jac)*wpg(ipg);
    end
    eleDofs = node2dof(Elementos(iele,:),nDofNod);
    R(eleDofs) = R(eleDofs) + Re;  
end

%% Reduccion Matriz
isFixed = reshape(bc',[],1);
isFree = ~isFixed;

Rr = R; 
% Rr = reshape(R',[],1);

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
stress = zeros(nel,nNodEle,4);
strain = zeros(nel,nNodEle,4);
uNod = [ -1 -1
          1 -1
          1  1
         -1  1 ];

Ni = shapefuns   (uNod,eleType);
dN = shapefunsder(uNod,eleType);     
for iele = 1:nel
    nodesEle = Nodos(Elementos(iele,:),:);
    for inode = 1:nNodEle
        % Derivadas de x,y, respecto de ksi, eta
        jac = dN(:,:,inode)*nodesEle;                      
        % Derivadas de las funciones de forma respecto de x,y.
        dNxy = jac\dN(:,:,inode);          % dNxy = inv(jac)*dN(:,:,ipg)
        r = Ni(:,:,inode)*nodesEle(:,1);
        
        B = zeros(size(C,2),nDofNod*nNodEle);
        B(1,1:2:7) = dNxy(1,:);
        B(2,1:2:7) = Ni(:,:,inode)/r;
        B(3,2:2:8) = dNxy(2,:);
        B(4,1:2:7) = dNxy(2,:);
        B(4,2:2:8) = dNxy(1,:);
        
        eleDofs = node2dof(Elementos(iele,:),nDofNod);
        def = B*D(eleDofs);
        strain(iele,inode,:) = def;
        stress(iele,inode,:) = C*def;
    end
end

%% Configuración deformada
D = (reshape(D,nDofNod,[]))';
nodePosition = Nodos + D(:,1:2);

%% Graficación
bandplot(Elementos,nodePosition,stress(:,:,2),[],'k',5);
meshplot(Elementos,Nodos,'b',0)

%% Solución analítica
a = Ro;
b = R1;
r = 0*Ro + 1*R1;

sr = (3 + NU)/8 * ( a^2 + b^2 - r^2 - (a*b/r)^2 ) * RHO * w^2;
st = (3 + NU)/8 * ( a^2 + b^2 - (1 + 3*NU)/(3 + NU)*r^2 + (a*b/r)^2 ) * RHO * w^2;
u  = (3 + NU)*(1 - NU)/(8*E) * ( a^2 + b^2 - (1 + NU)/(3 + NU)*r^2 + (1 + NU)/(1 - NU)*(a*b/r)^2 ) * RHO * w^2 * r;
ez = -NU/E * (sr + st);