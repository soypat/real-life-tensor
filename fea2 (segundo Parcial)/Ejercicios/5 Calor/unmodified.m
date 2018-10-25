clear
clc
close all

TestType = 'Sigma_y' ; % 'Sigma_y' ; % 'Sigma_x' ; %
eleType = 'Q4';% 'Q4', 'Q8', 'Q9'
escala=1; %Escala desplazamientos
%% Discretizacion
aux=load('elementosCQ4.txt');
elementos=aux(:,2:5);
aux=load('nodosCQ4.txt');
nodos=aux(:,2:3);

%% Propiedades del Material
E=200e9;
nu=0.3;
%Plane Strain
%C = (E/((1+nu)*(1-2*nu)))*[1-nu    nu      0;
                           %nu  1-nu      0;
                            %0    0  0.5-nu];
%Plane Stress
C = (E/(1-nu^2))*[1 nu 0;
                  nu 1 0;
                  0 0 (1-nu)/2];

%para un caso termico
%lo que se deforma es independiente del material??
k=25;%conduccion del material!
Ct=[k 0;0 k];
alfa=1.25e-5;
%% Definiciones
nDofNod = 2;                    % grados de libertad por nodo
nel = size(elementos,1);         % elementos
nNod = size(nodos,1);           % nodos
nNodEle = size(elementos,2);     % nodos por elemento
nDofTot = nDofNod*nNod;         % grados de libertad
nDims = size(nodos,2);          % dimensiones del problema
%si tengo un problema termico tengo que tener en cuenta algunas
%consideraciones, voy a tener una dimension menos en el Kt!!
nDofNodt=1;
nDofTott = nDofNodt*nNod;
%% Condiciones de borde
bc = false(nNod,nDofNod);       % Matriz de condiciones de borde
bc([1,5,11,15],1:2) = true;


%%
figure(1)
Meshplot(elementos,nodos,bc,'k',1)

%% Puntos de Gauss
rsInt = 3*ones(1,2);
[wpg, upg, npg] = gauss(rsInt);

%% Matriz de rigidez
Kt = zeros(nDofTott);
K=zeros(nDofTot);
A = 0;
jmin = 1E10;
for iele = 1:nel
    Ket = zeros(nDofNodt*nNodEle);
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
        %esto no lo uso ahora, pero si tengo un problema mecanico lo dejo!
        B = zeros(size(C,2),nDofNod*nNodEle);
        B(1,1:2:nDofNod*nNodEle-1) = dNxy(1,:);
        B(2,2:2:nDofNod*nNodEle) = dNxy(2,:);
        B(3,1:2:nDofNod*nNodEle-1) = dNxy(2,:);
        B(3,2:2:nDofNod*nNodEle) = dNxy(1,:);
        Bt=dNxy;%el B para cargas termicas
        Djac = det(jac);
        Ke = Ke + B'*C*B*wpg(ipg)*Djac;
        Ket = Ket + Bt'*Ct*Bt*wpg(ipg)*Djac;
        A = A + wpg(ipg)*Djac;
        if Djac < jmin
            jmin = Djac;
        end
    end
    eleDofs = node2dof(elementos(iele,:),nDofNod);
    eleDofst=node2dof(elementos(iele,:),1);%ahora lo hago para carga termica
    Kt(eleDofst,eleDofst) = Kt(eleDofst,eleDofst) + Ket;
    K(eleDofs,eleDofs) = K(eleDofs,eleDofs) + Ke;
end
%%
%%Distribucion de temperaturas, con conocidos y desc
Co=[1,6,5,10,11,12,13,14,15];
X=[2,3,4,7,8,9];
Kxx=Kt(X,X);Kxc=Kt(X,Co);
Tc=zeros(length(Co),1);
Tc([1,2,5])=100;%los nodos con respecto a Co
Rt=Kxc*Tc;
Tx=Kxx\-Rt;
T=zeros(nNod,1);
T(X)=Tx;
T(Co)=Tc;
%grafico la temp
for i=1:nel
    temp(i,:)=T(elementos(i,:));
end
figure(2)
bandplot(elementos,nodos,temp);
%%
%%cargas Termicas, copio lo de la K y cambio algunas cosas
R=zeros(nNod,nDofNod);
for iele = 1:nel
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
        %esto no lo uso ahora, pero si tengo un problema mecanico lo dejo!
        B = zeros(size(C,2),nDofNod*nNodEle);
        B(1,1:2:nDofNod*nNodEle-1) = dNxy(1,:);
        B(2,2:2:nDofNod*nNodEle) = dNxy(2,:);
        B(3,1:2:nDofNod*nNodEle-1) = dNxy(2,:);
        B(3,2:2:nDofNod*nNodEle) = dNxy(1,:);
        Djac = det(jac);
        R(elementos(iele,:),[1 2]) = R(elementos(iele,:),[1 2]) + reshape(-B'*C*(N*T(elementos(iele,:))*[1 1 0]')*alfa*wpg(ipg)*Djac,2,nNodEle)';
        A = A + wpg(ipg)*Djac;
        if Djac < jmin
            jmin = Djac;
        end
    end
end

% Reduccion Matriz
isFixed = reshape(bc',[],1);
isFree = ~isFixed;
Rr = reshape(R',[],1);
% Ktt=zeros(30);
% Kt(1:2:end,1:2:end)=Kt;
% Ktt(2:2:end,2:2:end)=Kt;
% Rr=zeros(nNod*nDofNod,2);

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
%%
%Flujo
%partimos como calcular tensiones
%%
%Tensiones y flujo
q=zeros(nNodEle*nDofNod,1);
stress = zeros(nNodEle,3,nel);
for iele = 1:nel
    nodesEle = nodos(elementos(iele,:),:);
    eleDofst=node2dof(elementos(iele,:),1);
    TnodosEle=T(elementos(iele,:));
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
        Bt=dNxy;
        eleDofs = node2dof(elementos(iele,:),nDofNod);
        stress(inode,:,iele) = C*(B*D(eleDofs)-(alfa*TnodosEle(inode)*[1 1 0]'));
        q([2*eleDofst(inode)-1,2*eleDofst(inode)])=-k*Bt*T(eleDofst);
    end
end

%%
%grafico flujo
%separo en x e y para graficar!!
qx=q(1:2:end);
qy=q(2:2:end);

figure(3)
%Meshplot(elementos,nodos,bc,'k',1)
quiver(nodos(:,1),nodos(:,2),qx,qy)
%% Configuración deformada
nodePosition = nodos + escala*(reshape(D,nDofNod,[]))';
figure(1)
Meshplot(elementos,nodePosition,bc,'r',0)
%% Gráfico
S=1;
figure(4)
bandplot(elementos,nodePosition,squeeze(stress(:,S,:))',[],'k');