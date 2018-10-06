clear
% clc
close all

TestType = 'Sigma_y' ; % 'Sigma_y' ; % 'Sigma_x' ; %
eleType = 'Q8';% 'Q4', 'Q9'
escala=1; %Escala desplazamientos
% Discretizacion
aux=load('Q8elements.txt');
elementos=aux(:,2:9);%porque es Q4
aux=load('Q8nodes.txt');
nodos=aux(:,2:3);
% elementos=[1 5 9 8;
%     5 2 6 9;
%     9 6 3 7;
%     8 9 7 4];
% 
% nodos = [ 0.00      0.00
%     1.00      0.00
%     1.00      1.00
%     0.00      1.00
%     0.50      0.00
%     1.00      0.50
%     0.50      1.00
%     0.00      0.50
%     0.50      0.50];
%% Propiedades del Material
E=1;
nu=0.3;
%Plane Strain
%C = (E/((1+nu)*(1-2*nu)))*[1-nu    nu      0;
                           %nu  1-nu      0;
                            %0    0  0.5-nu];
%Plane Stress, porque en este caso esta cargado en el plano medio, tengo
%una chapa en 2D
C = (E/(1-nu^2))*[1 nu 0;
                  nu 1 0;
                  0 0 (1-nu)/2];


%% Definiciones
%reescribimos la matriz elementos y nodos, sino lo cargo de adina
          
nDofNod = 2;                    % grados de libertad por nodo
nel = size(elementos,1);         % elementos
nNod = size(nodos,1);           % nodos
nNodEle = size(elementos,2);     % nodos por elemento
nDofTot = nDofNod*nNod;         % grados de libertad
nDims = size(nodos,2);          % dimensiones del problema

%% Condiciones de borde
bc = false(nNod,nDofNod);       % Matriz de condiciones de borde
bc(1,1:2) = true;
bc(9,2) = true;
%%
figure(1)
Meshplot(elementos,nodos,bc,'k',1)

%% Cargas


%% Puntos de Gauss
rsInt = 3*ones(1,2);
[wpg, upg, npg] = gauss(rsInt);

%% Matriz de rigidez
K = zeros(nDofTot);
A = 0;
jmin = 1E10;
int = zeros(nDofTot);
for iele = 1:nel
    
    Ke = zeros(nDofNod*nNodEle);
    nodesEle = nodos(elementos(iele,:),:);
    eleDofs = node2dof(elementos(iele,:),nDofNod);
    eleDofs = reshape(eleDofs',[],1);
    
    for ipg = 1:npg
        % Punto de Gauss
        ksi = upg(ipg,1);
        eta = upg(ipg,2);
        
        % Funciones de forma respecto de ksi, eta
        N = shapefuns([ksi eta],eleType);
        
%         Nm(1,1:2:7)=N;
%         Nm(2,2:2:8)=N;
        
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
        %hago la integral me va a servir despues para calcular las fuerzas
%         int(eleDofs,eleDofs)=int(eleDofs,eleDofs)+Nm'*Nm*det(jac);
        if Djac < jmin
            jmin = Djac;
        end
    end
    
    K(eleDofs,eleDofs) = K(eleDofs,eleDofs) + Ke;
end
%% Cargas

R = zeros(nNod,nDofNod);

% %Fuerzas de volumen
% fv=reshape([zeros(1,nNod);-9.81*2000*ones(1,nNod)],nDofTot,1);
% Fv=int*fv; Fv=reshape(Fv,nDofNod,nNod)';
% R=R+Fv;

% Fuerzas en superficies
q1 =@(x) 1/2;
% q2 = @(x) -550000*(1-x/52.5);
wallnodes=[19 14;
            20 16;
            3 10;
            7 13;
            1 9];
=[nodos(wallnodes(:,1),2) nodos(wallnodes(:,2),2)];
Lwall=abs(diff(aux));
Nwall=size(Lwall,1);
Nq8=[4 2 -1;2 16 2;-1 2 4];
Fv=@(qv,L) L/15*Nq8*qv;

for iele=1:Nwall
    pos=0;
    qv=[q(pos); q(pos); q(pos)];
    if iele==1
        Ll=Lwall(iele,1);
        Lr=Lwall(iele,2);
        R(wallnodes(iele,1),1)=-Ll/2*q1;
        R(wallnodes(iele,2),1)=Lr/2*q1;
        R(wallnodes(iele+1,1),1)=-Ll/2*q1;
        R(wallnodes(iele+1,2),1)=Lr/2*q1;
    elseif iele~=Nwall
        Ll=Lwall(iele,1);
        Lr=Lwall(iele,2);

    else
        Ll=Lwall(iele,1);
        Lr=Lwall(iele,2);

    end
        
end
for iele=1:nel
    nodesNum = elementos(iele,:);
    nodesEle = nodos(nodesNum,:);
    %cargas derecha
    if  max(nodesEle(:,1)) == 1
        L = abs(nodesEle(3,2)-nodesEle(2,2));
        Flat=L/6*[2 1;1 2]*[q1 q1]';
        R(nodesNum([2 3]),1)=R(nodesNum([2 3]),1)+Flat;
% %     cargas izquierda
    elseif min(nodesEle(:,1)) == 0 
        L = abs(nodesEle(4,2)-nodesEle(1,2));
        Flat=L/6*[2 1;1 2]*[-q1 -q1]';
        R(nodesNum([1 4]),1)=R(nodesNum([1 4]),1)+Flat;
    end
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
    
% Configuración deformada
nodePosition = nodos + escala*(reshape(D,nDofNod,[]))';
figure(1)
Meshplot(elementos,nodePosition,bc,'r',0)
% Gráfico
figure(2)
bandplot(elementos,nodePosition,stress(:,:,S)',[],'k');

