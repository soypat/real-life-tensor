clear
close all
format short g
clc

% Discretizacion
nodes = [ 0.0  0.0
          1.0  0.0
          2.0  0.0
          3.0  0.0
          0.0  1.5
          1.25 0.75
          2.0  1.0
          3.0  1.25
          0.0  3.5
          1.0  3.0
          2.5  2.5
          3.0  3.0
          0.0  5.0
          1.5  5.0
          2.0  5.0
          3.0  5.0];

elements = [1  2  6  5
            2  3  7  6
            3  4  8  7
            5  6 10  9
            6  7 11 10
            7  8 12 11
            9 10 14 13
           10 11 15 14
           11 12 16 15];

nDofNod = 2;                    % Número de grados de libertad por nodo
nNodEle = 4;                    % Número de nodos por elemento
nel = size(elements,1);         % Número de elementos
nNod = size(nodes,1);           % Número de nodos
nDofTot = nDofNod*nNod;         % Número de grados de libertad

bc = false(nNod,nDofNod);       % Matriz de condiciones de borde
bc(1,1:2) = true;
bc(2:4,2) = true;

R = zeros(nNod,nDofNod);        % Vector de cargas
R(13,1)=(nodes(14,1)/2-nodes(13,1)/2);
R(14,1)=(nodes(14,1)/2-nodes(13,1)/2)+(nodes(15,1)/2-nodes(14,1)/2);
R(15,1)=(nodes(15,1)/2-nodes(14,1)/2)+(nodes(16,1)/2-nodes(15,1)/2);
R(16,1)=(nodes(16,1)/2-nodes(15,1)/2);
R(1,1)=-(nodes(2,1)/2-nodes(1,1)/2);
R(2,1)=-(nodes(2,1)/2-nodes(1,1)/2)-(nodes(3,1)/2-nodes(2,1)/2);
R(3,1)=-(nodes(3,1)/2-nodes(2,1)/2)-(nodes(4,1)/2-nodes(3,1)/2);
R(4,1)=-(nodes(4,1)/2-nodes(3,1)/2);
R(1,2)=-(nodes(5,2)/2-nodes(1,2)/2);
R(5,2)=-(nodes(5,2)/2-nodes(1,2)/2)-(nodes(9,2)/2-nodes(5,2)/2);
R(9,2)=-(nodes(9,2)/2-nodes(5,2)/2)-(nodes(13,2)/2-nodes(9,2)/2);
R(13,2)=-(nodes(13,2)/2-nodes(9,2)/2)+(nodes(14,1)-nodes(13,1));
R(4,2)=(nodes(8,2)/2-nodes(4,2)/2);
R(8,2)=(nodes(8,2)/2-nodes(4,2)/2)+(nodes(12,2)/2-nodes(8,2)/2);
R(12,2)=(nodes(12,2)/2-nodes(8,2)/2)+(nodes(16,2)/2-nodes(12,2)/2);
R(16,2)=(nodes(16,2)/2-nodes(12,2)/2)+(nodes(16,1)-nodes(15,1));
R(14,2)=(nodes(14,1)-nodes(13,1))+(nodes(15,1)-nodes(14,1));
R(15,2)=(nodes(15,1)-nodes(14,1))+(nodes(16,1)-nodes(15,1));


%% Matriz Constitutiva (plane stress)
% Propiedades del material
E = 10;
NU = -0.1;
C = E/(1 - NU^2)*[ 1.0     NU         0.0
                    NU    1.0         0.0
                   0.0    0.0     (1 - NU)/2 ];

%% Matriz de rigidez
% Gauss
[wpg, upg, npg] = gauss([2 2]);
upg = upg([1 3 4 2],:);
K = zeros(nDofTot);
nodeDofs = reshape(1:nDofTot,nDofNod,nNod)';
detjac = zeros(nel,nNodEle);
relArea = zeros(nel,1);
for iele = 1:nel
    Ke = zeros(nDofNod*nNodEle);
    nodesEle = nodes(elements(iele,:),:);
    for ipg = 1:npg
        % Punto de Gauss
        ksi = upg(ipg,1);
        eta = upg(ipg,2);  
        % Derivadas de las funciones de forma respecto de ksi, eta
        dN = shapefunsder([ksi eta],'Q4');  
        % Derivadas de x,y, respecto de ksi, eta
        jac = dN*nodesEle;                      
        % Derivadas de las funciones de forma respecto de x,y.
        dNxy = jac\dN;          % dNxy = inv(jac)*dN
        
        B = zeros(size(C,2),nDofNod*nNodEle);
        B(1,1:2:nNodEle*2-1) = dNxy(1,:);
        B(2,2:2:nNodEle*2)   = dNxy(2,:);
        B(3,1:2:nNodEle*2-1) = dNxy(2,:);
        B(3,2:2:nNodEle*2)   = dNxy(1,:); 

        Ke = Ke + B'*C*B*wpg(ipg)*det(jac);
        
        % Alacenado de jacobiano
        detjac(iele,ipg) = det(jac);
        relArea(iele) = relArea(iele) + wpg(ipg) * det(jac);
    end
    eleDofs = nodeDofs(elements(iele,:),:);
    eleDofs = reshape(eleDofs',[],1);
    K(eleDofs,eleDofs) = K(eleDofs,eleDofs) + Ke;  
end

dist = min(detjac,[],2).*relArea;
[maxdist, maxdistele] = max(dist)

%% Reduccion Matriz
isFixed = reshape(bc',[],1);
isFree = ~isFixed;

Rr = reshape(R',[],1);

% Solver
Dr = K(isFree,isFree)\Rr(isFree);

% Reconstruccion
D = zeros(nDofTot,1);
D(isFree) = D(isFree) + Dr;

% Reacciones
Rv = K(isFixed,isFree)*D(isFree);
reacciones = nan(nDofTot,1);
reacciones(isFixed) = Rv;
reacciones = (reshape(reacciones,nDofNod,[]))';

%% Graficos
figure
meshplot(elements,nodes,'b')
hold on
Dreshaped = (reshape(D,nDofNod,[]))';
nodePosition = nodes + Dreshaped(:,1:2);
meshplot(elements,nodePosition,'r')

%% Recuperacion de tensiones en los nodos
stressNod = zeros(nel,nNodEle,3);
uNod = [-1 -1
         1 -1
         1  1
        -1  1];
for iele = 1:nel
    nodesEle = nodes(elements(iele,:),:);
    for inode = 1:nNodEle
        % Punto de Gauss
        ksi = uNod(inode,1);
        eta = uNod(inode,2);  
        % Derivadas de las funciones de forma respecto de ksi, eta
        dN = shapefunsder([ksi eta],'Q4');  
        % Derivadas de x,y, respecto de ksi, eta
        jac = dN*nodesEle;                      
        % Derivadas de las funciones de forma respecto de x,y.
        dNxy = jac\dN;          % dNxy = inv(jac)*dN
        
        B = zeros(size(C,2),nDofNod*nNodEle);
        B(1,1:2:nNodEle*2-1) = dNxy(1,:);
        B(2,2:2:nNodEle*2)   = dNxy(2,:);
        B(3,1:2:nNodEle*2-1) = dNxy(2,:);
        B(3,2:2:nNodEle*2)   = dNxy(1,:); 
        
        eleDofs = nodeDofs(elements(iele,:),:);
        eleDofs = reshape(eleDofs',[],1);
        stressNod(iele,inode,:) = C*B*D(eleDofs);
    end
end

figure
bandplot(elements,nodes,stressNod(:,:,3),[],'k')

%% Tensiones en los puntos de superconvergencia
stressSuper = zeros(nel,nNodEle,3);
for iele = 1:nel
    nodesEle = nodes(elements(iele,:),:);
    for inode = 1:nNodEle
        % Punto de Gauss
        ksi = 0;
        eta = 0;  
        % Derivadas de las funciones de forma respecto de ksi, eta
        dN = shapefunsder([ksi eta],'Q4'); 
        % Derivadas de x,y, respecto de ksi, eta
        jac = dN*nodesEle;                      
        % Derivadas de las funciones de forma respecto de x,y.
        dNxy = jac\dN;          % dNxy = inv(jac)*dN
        
        B = zeros(size(C,2),nDofNod*nNodEle);
        B(1,1:2:nNodEle*2-1) = dNxy(1,:);
        B(2,2:2:nNodEle*2)   = dNxy(2,:);
        B(3,1:2:nNodEle*2-1) = dNxy(2,:);
        B(3,2:2:nNodEle*2)   = dNxy(1,:);
        
        eleDofs = nodeDofs(elements(iele,:),:);
        eleDofs = reshape(eleDofs',[],1);
        stressSuper(iele,inode,:) = C*B*D(eleDofs);
    end
end

%% Promediado de tensiones en los nodos
avgStress = zeros(nNod,3);
for inode = 1:nNod
    [I,J] = find(elements == inode);
    nShare = length(I);
    for ishare = 1:nShare
        avgStress(inode,:) = avgStress(inode,:) + squeeze(stressSuper(I(ishare),J(ishare),:))';
    end
    avgStress(inode,:) = avgStress(inode,:) / nShare;
end

%% calculo de eta_el, e2, U2
invC = C\eye(3);
eta_el = zeros(nel,1);
e2_el = zeros(nel,1);
U2_el = zeros(nel,1);
for iele = 1:nel
    nodesEle = nodes(elements(iele,:),:);
    for inode = 1:nNodEle
        for ipg = 1:npg
            % Punto de Gauss
            ksi = upg(ipg,1);
            eta = upg(ipg,2);

            % Derivadas de las funciones de forma respecto de ksi, eta
            dN = shapefunsder([ksi eta],'Q4');
            % Derivadas de x,y, respecto de ksi, eta
            jac = dN*nodesEle;                      
            % Derivadas de las funciones de forma respecto de x,y.
            dNxy = jac\dN;          % dNxy = inv(jac)*dN

            % funciones de forma
            N = shapefuns([ksi eta],'Q4');
            eleStress =  squeeze(stressNod(iele,inode,:));    % tensiones "directas"
            starStress = ( N * avgStress(elements(iele,:),:) )';    % tensiones mejoradas
            
            e2_el(iele) = e2_el(iele) + (starStress - eleStress)' * ... 
                        invC * (starStress - eleStress) * wpg(ipg) * det(jac);

            U2_el(iele) = U2_el(iele) + eleStress' * invC * eleStress * ...
                          wpg(ipg) * det(jac);
        end
    end
    eta_el(iele) = sqrt( e2_el(iele) / (e2_el(iele) + U2_el(iele)) );
    
end

etaG = sqrt( sum(e2_el) / (sum(e2_el) + sum(U2_el)) );

disp('Error global por energía de deformación: ')
disp (etaG)
fprintf('\n')
disp('Error elemental por energía de deformación: ')
disp(([1:nel; eta_el'])')

%% Pregunta 2
K2 = K;
R2 = reshape(R',[],1);
bc2 = reshape(bc',[],1);
[i,j]=find((abs(nodes(:,1)-3) <= 1e-3) & (abs(nodes(:,2)-0) <= 1e-3));
masterDof = nodeDofs(i,j);
c = 0;
for iele = 1:nel
    nodesEle = nodes(elements(iele,:),:);
    if abs(nodesEle(2,1)-3) <= 1e-3
        retained = zeros(size(K2,1),1);
        eleDofs = nodeDofs(elements(iele,3),1);
        eleDofs = reshape(eleDofs',[],1); % Si fuese mas de un valor sirve esto
        retained([eleDofs masterDof]) = [1; -1];
        K2 = [K2 retained
              retained' 0];
        R2 = [R2; 0];
        bc2 = [bc2; 0];
        c = c + 1;
    end
end

%% Reduccion Matriz

isFixed = logical(bc2);
isFree = ~isFixed;

% Solver
Dr2 = K2(isFree,isFree)\R2(isFree);

% Reconstruccion
D2 = zeros(nDofTot+c,1);
D2(isFree) = D2(isFree) + Dr2;

% Reacciones
Rv2 = K2(isFixed,isFree)*D2(isFree);
reacciones2 = nan(nDofTot+c,1);
reacciones2(isFixed) = Rv2;
reacciones2 = (reshape(reacciones2(1:end-c),nDofNod,[]))';

%% Graficos
figure
meshplot(elements,nodes,'b')
hold on
Dreshaped2 = (reshape(D2(1:end-c),nDofNod,[]))';
nodePosition2 = nodes + Dreshaped2(:,1:2);
meshplot(elements,nodePosition2,'r')

%% Determinación de la fuerza
F = [];
for iele = 1:nel
    nodesEle = nodes(elements(iele,:),:);
    if abs(nodesEle(2,1)-3) <= 1e-3
        eleDofs = nodeDofs(elements(iele,3),1);
        eleDofs = reshape(eleDofs',[],1); % Si fuese mas de un valor sirve esto
        F = [F D2(eleDofs)];
    end
end

RR = K*D2(1:end-c)








