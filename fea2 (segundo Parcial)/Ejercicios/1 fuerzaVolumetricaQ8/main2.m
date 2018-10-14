%{
Preguntas para Sebas: Como calcular fuerzas volumetricas. Vale usar 4|J|=A?
Que nos dan para el parcial?

%}
escala=1;
aux=load('elementosEjemplo.txt');
elementos=aux(:,2:9);%porque es Q8
aux=load('nodEjemplo.txt');
nodos=aux(:,2:3); %mm
numeracion=aux(:,1);
[nodos, elementos]=nodekill(nodos,numeracion,elementos);
% nodxs=-1*ones(max(nodenumbering),2); %Creo una nueva matriz que tiene los nodos en la fila correspondiente a su numero asignado
% for i=1:size(nodos,1)
%     nodxs(nodenumbering(i),[1 2])=nodos(i,:); %Tengo que hacer esto porque Mati o quien sea programo esta porqueria para que tome los nodos segun su linea. Parece a proposito este engendro de paradigma de programacion, seguro lo es.
% end
% fakeNodeNumber=max(nodenumbering);
% fNN=fakeNodeNumber;
% nodos=nodxs; %Finalmente asigno la matriz nodos
% Seguro podria hacer el programa de cero, pero para que? No quiero
% aprender.
ticStart=tic;
%% Datellis
E=10; %N/mm^2
nu=0.3;
C = (E/((1+nu)*(1-2*nu)))*[1-nu    nu      0;
                           nu  1-nu      0;
                            0    0  0.5-nu];
%Plain Stress
          % tal cual como me lo dieron a mi:

nDofNod = 2;                    % grados de libertad por nodo
nel = size(elementos,1);         % elementos
nNod = size(nodos,1);           % nodos
nNodEle = size(elementos,2);     % nodos por elemento
nDofTot = nDofNod*nNod;         % grados de libertad
nDims = size(nodos,2);          % dimensiones del problema
%% CB
bc = false(nNod,nDofNod);       % Matriz de condiciones de borde
bc([1 3],:)=true;
%% PLOT
figure(1)
myMeshplot(elementos,nodos,bc,'k',1,1) %Eligo si quiero con numeración con el ultimo parametro (1/0)

%% Puntos de Gauss
rsInt = 3*ones(1,2);
[wpg, upg, npg] = gauss(rsInt);
Nm=zeros(2,16);
%% Matriz de rigidez
K = zeros(nDofTot);
KT = zeros(nDofTot/2);

A = 0;
jmin = 1E10;
inti = zeros(nDofTot);
Areas=zeros(nel,1);

for iele = 1:nel
    nodesEle = nodos(elementos(iele,:),:);
    %Mecánicas
    Ke = zeros(nDofNod*nNodEle);
    eleDofs = node2dof(elementos(iele,:),nDofNod);
    eleDofs = reshape(eleDofs',[],1);
    %Térmicas
%     KTe= zeros(1*nNodEle);
%     eleDofsT=node2dof(elementos(iele,:),1);
%     eleDofsT=reshape(eleDofsT',[],1);
    for ipg = 1:npg
        % Punto de Gauss
        ksi = upg(ipg,1);
        eta = upg(ipg,2);
        
        % Funciones de forma respecto de ksi, eta
        N = shapefuns([ksi eta],'Q8');
        
        Nm(1,1:2:15)=N;
        Nm(2,2:2:16)=N;
        
        % Derivadas de las funciones de forma resp16to de ksi, eta
        dN = shapefunsder([ksi eta],'Q8');
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
        %Mecánica
        Ke = Ke + B'*C*B*wpg(ipg)*Djac;
        %Térmica
%         KTe = KTe + BT'*Ct*BT*wpg(ipg)*Djac;
        
        A = A + wpg(ipg)*Djac;

        if Djac < jmin
            jmin = Djac;
        end
    end
%     KT(eleDofsT,eleDofsT)=KT(eleDofsT,eleDofsT)+KTe;
    Areas(iele)=A;
    K(eleDofs,eleDofs) = K(eleDofs,eleDofs) + Ke;
    A=0;
end
%% Cargas
R = zeros(nNod,nDofNod);
g=1; %m/s
t=1; %espesor
rho=1;
f=rho*g;
for iele=1:nel
    nodesEle = nodos(elementos(iele,:),:);
    index=elementos(iele,:);
    for ipg=1:npg
        ksi = upg(ipg,1);
        eta = upg(ipg,2);
        
        weight=wpg(ipg);
        
        N = shapefuns([ksi eta],'Q8');
        dN = shapefunsder([ksi eta],'Q8');
        
        
        jac = dN*nodesEle;
        
        dNxy = jac\dN;          % dNxy = inv(jac)*dN

        B = zeros(size(C,2),nDofNod*nNodEle);
        B(1,1:2:nDofNod*nNodEle-1) = dNxy(1,:);
        B(2,2:2:nDofNod*nNodEle) = dNxy(2,:);
        B(3,1:2:nDofNod*nNodEle-1) = dNxy(2,:);
        B(3,2:2:nDofNod*nNodEle) = dNxy(1,:);
        
        Detjac=det(jac);
        
        R(index,2)=R(index,2)+N'*f*Detjac*weight;
    end
end

%% Reducción y Resolución
isFixed = reshape(bc',[],1);
isFree = ~isFixed;

Rr = reshape(R',[],1);

% Solver
ticSolve=tic;
Dr = K(isFree,isFree)\Rr(isFree);
tocSolve=toc(ticSolve);
% Reconstrucción
D = zeros(nDofTot,1);
D(isFree) = D(isFree) + Dr;

% Reacciones
Rv = K(isFixed,isFree)*D(isFree);
reacciones = nan(nDofTot,1);
reacciones(isFixed) = Rv;
reacciones = (reshape(reacciones,nDofNod,[]))';

%% Tensiones y graficos
 uNod = [-1 -1
             1 -1
             1  1
            -1  1
             0 -1
             1  0
             0  1
            -1  0
             0  0];
stress = zeros(nNodEle,nel,3);
for iele = 1:nel
    nodesEle = nodos(elementos(iele,:),:);
    for inode = 1:nNodEle
        % Punto de Gauss
        ksi = uNod(inode,1);
        eta = uNod(inode,2);
        % Derivadas de las funciones de forma respecto de ksi, eta
        dN  = shapefunsder([ksi eta],'Q8');
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
sigx=stress(:,:,1);
sigy=stress(:,:,2);
sigxy=stress(:,:,3);
sigvm=sqrt(sigx.^2+sigy.^2 -sigx.*sigy + 3*sigxy.^2);
S=2;
    
% Configuración defo+rmada
nodePosition = nodos + escala*(reshape(D,nDofNod,[]))';
figure(1)
myMeshplot(elementos,nodePosition,bc,'r',0,1)
% Gráfico
figure(2)
bandplot(elementos,nodePosition,sigvm',[],'k');
tocEnd=toc(ticStart);
fprintf('El programa termino en %0.2f segundos\nLa solución fue encontrada en %0.3f segundos\n',tocEnd,tocSolve)