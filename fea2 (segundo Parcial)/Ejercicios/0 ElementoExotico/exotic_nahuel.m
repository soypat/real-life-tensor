clear
clc
close all

TestType = 'Sigma_y' ; % 'Sigma_y' ; % 'Sigma_x' ; %
eleType = 'Q5';% 'Q4', 'Q8', 'Q9'
escala=1; %Escala desplazamientos
%% Discretizacion
e=1;
Nn=5;
syms x y
b=1; h=1;
nodos=[[-b -h];[b -h];[b h];[-b h];[0 0]];
elementos=[1 2 3 4 5];
A=sym('A',[5 5]);
X=[1 x y x*y x^2*y^2];%armo el triangulo de tartaglia
for i=1:Nn
A(i,:)=subs(X,[x y],nodos(i,:));
end
N=X*inv(A);

for k=1:Nn
    Nu=N(k);
    syms x y
    Nx=diff(Nu,x);
    Ny=diff(Nu,y);
    B(1,2*k-1)=Nx;
    B(3,2*k-1)=Ny;
    B(2,2*k)=Ny;
    B(3,2*k)=Nx;    
end
 E=200e9; nu=0.3;
% C=(E/(1-nu^2))*[1 nu 0;nu 1 0;0 0 (1-nu)/2];
C = (E/((1+nu)*(1-2*nu)))*[1-nu    nu      0;
                           nu  1-nu      0;
                            0    0  0.5-nu];
BCB=B.'*C*B;
K=int(BCB,x,-b,b);
K=int(K,y,-h,h);
K=double(K);

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
R = zeros(nNod,nDofNod);        
% Vector de cargas
R(3,2)=2;
R(4,2)=-2;

%% Puntos de Gauss
rsInt = 3*ones(1,2);
[wpg, upg, npg] = gauss(rsInt);

%%
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
if strcmp(eleType,'Q5')
    uNod = [-1 -1
             1 -1
             1  1
            -1  1
             0  0];
elseif strcmp(eleType,'Q4')
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
stress = zeros(nNodEle,3,nel);
for iele = 1:nel
    nodesEle = nodos(elementos(iele,:),:);
    for inode = 1:nNodEle
        B=subs(B,{x,y},{nodesEle(inode,1),nodesEle(inode,1)});
        B=double(B);
        eleDofs = node2dof(elementos(iele,:),nDofNod);
        stress(inode,:,iele) = C*B*D(eleDofs);
    end
end
Svm=sqrt(stress(:,1).^2-stress(:,1).*stress(:,2)+stress(:,2).^2+3*stress(:,3).^2);
% S=1;
% elementos=[1 2 3 4];    
% % Configuración deformada
% nodePosition = nodos + escala*(reshape(D,nDofNod,[]))';
% figure(1)
% Meshplot(elementos,nodePosition,bc,'r',0)
% % Gráfico
% figure(2)
% bandplot(elementos,nodePosition,stress(:,:,S)',[],'k');

