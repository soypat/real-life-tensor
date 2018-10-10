t=1;
escala=1;
aux=load('elementos2.txt');
elementos=aux(:,2:9);%porque es Q8
aux=load('nodos2.txt');
nodos=aux(:,2:3); %mm
numeracion=aux(:,1);
[nodos, elementos]=nodekill(nodos,numeracion,elementos);
ticStart=tic;
%% Datellis
E=200e9; %N/mm^2
nu=0.3;
C = (E/(1-nu^2))*[1 nu 0;
                  nu 1 0;
                  0 0 (1-nu)/2]; %Plain Stress
          % tal cual como me lo dieron a mi:
Ct=17*[1 0;0 1]; %W/mK
alfa=1.2e-6;
%Termico

nDofNod = 2;                    % grados de libertad por nodo
nel = size(elementos,1);         % elementos
nNod = size(nodos,1);           % nodos
nNodEle = size(elementos,2);     % nodos por elemento
nDofTot = nDofNod*nNod;         % grados de libertad
nDims = size(nodos,2);          % dimensiones del problema
%% CB
bc = false(nNod,nDofNod);       % Matriz de condiciones de borde
for i=1:nNod
    if nodos(i,2)==0
        bc(i,:)=true;
    end
end
%% PLOT
figure(1)
myMeshplot(elementos,nodos,bc,'k',1,1) %Eligo si quiero con numeraci�n con el ultimo parametro (1/0)

%% Puntos de Gauss
rsInt = 3*ones(1,2);
[wpg, upg, npg] = gauss(rsInt);

%% Matriz de rigidez
K = zeros(nDofTot);
KT = zeros(nDofTot/2);
A = 0;
jmin = 1E10;
int = zeros(nDofTot);
Areas=zeros(nel,1);
for iele = 1:nel
    nodesEle = nodos(elementos(iele,:),:);
    %Mec�nicas
    Ke = zeros(nDofNod*nNodEle);
    eleDofs = node2dof(elementos(iele,:),nDofNod);
    eleDofs = reshape(eleDofs',[],1);
    %T�rmicas
%     KTe= zeros(1*nNodEle);
%     eleDofsT=node2dof(elementos(iele,:),1);
%     eleDofsT=reshape(eleDofsT',[],1);
    for ipg = 1:npg
        % Punto de Gauss
        ksi = upg(ipg,1);
        eta = upg(ipg,2);
        
        % Funciones de forma respecto de ksi, eta
        N = shapefuns([ksi eta],'Q8');
        
%         Nm(1,1:2:7)=N;
%         Nm(2,2:2:8)=N;
        
        % Derivadas de las funciones de forma respecto de ksi, eta
        dN = shapefunsder([ksi eta],'Q8');
        % Derivadas de x,y, respecto de ksi, eta
        jac = dN*nodesEle;
        % Derivadas de las funciones de forma respecto de x,y.
        dNxy = jac\dN;          % dNxy = inv(jac)*dN
%         %T�rmica
%         BT = zeros(2,4);
%         BT(1,:) = dNxy(1,:);
%         BT(2,:) = dNxy(2,:);
        %Mec�nica
        B = zeros(size(C,2),nDofNod*nNodEle);
        B(1,1:2:nDofNod*nNodEle-1) = dNxy(1,:);
        B(2,2:2:nDofNod*nNodEle) = dNxy(2,:);
        B(3,1:2:nDofNod*nNodEle-1) = dNxy(2,:);
        B(3,2:2:nDofNod*nNodEle) = dNxy(1,:);
        
        Djac = det(jac);
        Areas(iele)=abs(Djac*4);
        %Mec�nica
        Ke = Ke + B'*C*B*wpg(ipg)*Djac;
        %T�rmica
%         KTe = KTe + BT'*Ct*BT*wpg(ipg)*Djac;
        
        A = A + wpg(ipg)*Djac;
        %hago la integral me va a servir despues para calcular las fuerzas
%         int(eleDofs,eleDofs)=int(eleDofs,eleDofs)+Nm'*Nm*det(jac);
        if Djac < jmin
            jmin = Djac;
        end
    end
%     KT(eleDofsT,eleDofsT)=KT(eleDofsT,eleDofsT)+KTe;
    K(eleDofs,eleDofs) = K(eleDofs,eleDofs) + Ke;
end
%% Cargas
R = zeros(nNod,nDofNod);

%Fuerzas de volumen (Fuerza centrifuga!!!!)


%Carga Lineal
% Fuerzas en superficies
% q1 =@(x) -2*(150 - (210-x))/150; %MPa or N/mm^2
% q2 = @(x) -550000*(1-x/52.5);
% wallnodes=load('loadnod.txt');
wnodespos=nodos(wallnodes,2);
Lwall=diff(wnodespos);
Nwall=size(wallnodes,1);
Nq8=[4 2 -1;2 16 2;-1 2 4];
Fv=@(qv,L) L/15*Nq8*qv;

for iele=1:Nwall
    if iele<=(Nwall-1)/2 %para que no salga de indice
        index=[iele*2-1 iele*2 iele*2+1];
        pos=wnodespos(index);
        qv=[q1(pos(1)); q1(pos(2)) ;q1(pos(3))];
        rvl=(Lwall(iele*2,1)/15)*Nq8*qv;
        R(wallnodes(iele*2-1,1),1)=R(wallnodes(iele*2-1,1),1)+rvl(1);
        R(wallnodes(iele*2,1),1)=R(wallnodes(iele*2,1),1)+rvl(2);
        R(wallnodes(iele*2+1,1),1)=R(wallnodes(iele*2+1,1),1)+rvl(3);
        
    end
    %No esta definido para otro tipo de elemento que no sea Q8
end
%% Reducci�n y Resoluci�n
isFixed = reshape(bc',[],1);
isFree = ~isFixed;

Rr = reshape(R',[],1);

% Solver
ticSolve=tic;
Dr = K(isFree,isFree)\Rr(isFree);
tocSolve=toc(ticSolve);
% Reconstrucci�n
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
S=1;
    
% Configuraci�n deformada
nodePosition = nodos + escala*(reshape(D,nDofNod,[]))';
figure(1)
Meshplot(elementos,nodePosition,bc,'r',0)
% Gr�fico
figure(2)
bandplot(elementos,nodePosition,stress(:,:,S)',[],'k');
tocEnd=toc(ticStart);
fprintf('El programa termino en %0.2f segundos\nLa soluci�n fue encontrada en %0.3f segundos\n',tocEnd,tocSolve)