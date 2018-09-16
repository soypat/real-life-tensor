
clear
clc
close all
format short g

nodes = load('PA_Nod_Q8.txt');
nodes(:,[1 4 5]) = [];

elements =load('PA_Elem_Q8.txt');
elements(:,[1 10]) = [];
elements = elements(:,[3 4 1 2 7 8 5 6]);

nDofNod = 2;                    % Numero de grados de libertad por nodo
nNodEle = size(elements,2);     % Numero de nodos por elemento
nel = size(elements,1);         % Numero de elementos
nNod = size(nodes,1);           % Numero de nodos
nDofTot = nDofNod*nNod;         % Numero de grados de libertad

bc = false(nNod,nDofNod);       % Matriz de condiciones de borde
bc(nodes(:,1) == 0,1) = true;
bc(nodes(:,2) == 0,2) = true;

% meshplot(elements,nodes,'b')

%% Gauss           
a   = sqrt(.6);
% Ubicaciones puntos de Gauss
upg = [ -a  -a
         a  -a
         a   a
        -a   a 
         0  -a
         a   0
         0   a
        -a   0
         0   0];    
% Numero de puntos de Gauss
npg = size(upg,1);
wpg = [ones(1,4)*25/81 ones(1,4)*40/81 64/81];

%% Matriz Constitutiva (plane stress)
t = 1;
E = 1;
NU = 0.3;
C = E/(1 - NU^2)*[ 1.0     NU         0.0
                    NU    1.0         0.0
                   0.0    0.0     (1 - NU)/2 ];
               
%% Matriz de rigidez
K = zeros(nDofTot);
nodeDofs = reshape(1:nDofTot,nDofNod,nNod)';
for iele = 1:nel
    Ke = zeros(nDofNod*nNodEle);
    nodesEle = nodes(elements(iele,:),:);
    for ipg = 1:npg
        % Punto de Gauss
        ksi = upg(ipg,1);
        eta = upg(ipg,2);  
        % Derivadas de las funciones de forma respecto de ksi, eta
        dN = shapefunsder([ksi eta],'Q8');  
        % Derivadas de x,y, respecto de ksi, eta
        jac = dN*nodesEle;                      
        % Derivadas de las funciones de forma respecto de x,y.
        dNxy = jac\dN;          % dNxy = inv(jac)*dN
        
        B = zeros(size(C,2),nDofNod*nNodEle);
        B(1,1:2:15) = dNxy(1,:);
        B(2,2:2:16) = dNxy(2,:);
        B(3,1:2:15) = dNxy(2,:);
        B(3,2:2:16) = dNxy(1,:); 

        Ke = Ke + B'*C*B*wpg(ipg)*det(jac);
    end
    eleDofs = nodeDofs(elements(iele,:),:);
    eleDofs = reshape(eleDofs',[],1);
    K(eleDofs,eleDofs) = K(eleDofs,eleDofs) + Ke;  
end

%% Carga
F = 1e-2; % [N/m]
% Ubicaciones puntos de Gauss
upg = [-a 0 a];
% Numero de puntos de Gauss
npg = size(upg,2);
wpg = [5/9 8/9 5/9];
R = zeros(nNod,nDofNod);        % Vector de cargas
nodeDofs = reshape(1:nDofTot,nDofNod,nNod)';
for iele = 1:nel
    nodesEle = nodes(elements(iele,:),:);
    if nodesEle(3,2) >= 15-10e-3
        for ipg = 1:npg
            % Punto de Gauss
            ksi = upg(ipg);
            eta = 1;  
            % Derivadas de las funciones de forma respecto de ksi, eta
            dN = shapefunsder([ksi eta],'Q8'); 
            % Derivadas de x,y, respecto de ksi, eta
            jac = dN*nodesEle;                     
            % Derivadas de las funciones de forma respecto de x,y.
            dNxy = jac\dN;          % dNxy = inv(jac)*dN
            N = shapefuns([ksi eta],'Q8');
            N = N([4 7 3]);
            R(elements(iele,[4 7 3]),2) = R(elements(iele,[4 7 3]),2) + ...
                N'*N*F*ones(3,1)*jac(1,1)*wpg(ipg)*t;
        end
    end
end
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

%% Recuperacion de tensiones en los nodos
stress = zeros(nel,nNodEle,3);
uNod = [-1 -1
         1 -1
         1  1
        -1  1
         0 -1
         1  0
         0  1
        -1  0];
for iele = 1:nel
    nodesEle = nodes(elements(iele,:),:);
    for inode = 1:nNodEle
        % Punto de Gauss
        ksi = uNod(inode,1);
        eta = uNod(inode,2);  
        % Derivadas de las funciones de forma respecto de ksi, eta
        dN = shapefunsder([ksi eta],'Q8');  
        % Derivadas de x,y, respecto de ksi, eta
        jac = dN*nodesEle;                      
        % Derivadas de las funciones de forma respecto de x,y.
        dNxy = jac\dN;          % dNxy = inv(jac)*dN
        
        B = zeros(size(C,2),nDofNod*nNodEle);
        B(1,1:2:15) = dNxy(1,:);
        B(2,2:2:16) = dNxy(2,:);
        B(3,1:2:15) = dNxy(2,:);
        B(3,2:2:16) = dNxy(1,:);
        
        eleDofs = nodeDofs(elements(iele,:),:);
        eleDofs = reshape(eleDofs',[],1);
        stress(iele,inode,:) = C*B*D(eleDofs);
    end
end

%% tensiones en los puntos de superconvergencia
stressSuper = zeros(nel,4,3);
uNod = [-1  -1
         1  -1
         1   1
        -1   1]/sqrt(3);
for iele = 1:nel
    nodesEle = nodes(elements(iele,:),:);
    for inode = 1:4
        % Punto de Gauss
        ksi = uNod(inode,1);
        eta = uNod(inode,2);  
        % Derivadas de las funciones de forma respecto de ksi, eta
        dN = shapefunsder([ksi eta],'Q8');
        % Derivadas de x,y, respecto de ksi, eta
        jac = dN*nodesEle;                      
        % Derivadas de las funciones de forma respecto de x,y.
        dNxy = jac\dN;          % dNxy = inv(jac)*dN
        
        B = zeros(size(C,2),nDofNod*nNodEle);
        B(1,1:2:15) = dNxy(1,:);
        B(2,2:2:16) = dNxy(2,:);
        B(3,1:2:15) = dNxy(2,:);
        B(3,2:2:16) = dNxy(1,:);
        
        eleDofs = nodeDofs(elements(iele,:),:);
        eleDofs = reshape(eleDofs',[],1);
        stressSuper(iele,inode,:) = C*B*D(eleDofs);
    end
end

% Extrapolo a los Nodos
a = sqrt(3);
rsExt = a*[-1  -1
            1  -1
            1   1
           -1   1
            0  -1
            1   0
            0   1
           -1   0];
stressExtra = zeros(nel,nNodEle,3);
for iele = 1:nel
    for inode = 1:nNodEle
        r = rsExt(inode,1);
        s = rsExt(inode,2);
        
        N = shapefuns([r s],'Q4');
        
        stressExtra(iele,inode,:) = N * squeeze(stressSuper(iele,:,:));
    end
end

%% Gauss           
[wpg, upg, npg] = gauss([3 3]);

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
            dN = shapefunsder([ksi eta],'Q8');
            % Derivadas de x,y, respecto de ksi, eta
            jac = dN*nodesEle;                      
            % Derivadas de las funciones de forma respecto de x,y.
            dNxy = jac\dN;          % dNxy = inv(jac)*dN

            % funciones de forma
            N = shapefuns([ksi eta],'Q8');
            eleStress =  squeeze(stress(iele,inode,:));            % tensiones "directas"
            starStress = squeeze(stressExtra(iele,inode,:)); % tensiones mejoradas

            e2_el(iele) = e2_el(iele) + (starStress - eleStress)' * ... 
                        invC * (starStress - eleStress) * wpg(ipg) * det(jac);

            U2_el(iele) = U2_el(iele) + eleStress' * invC * eleStress * ...
                          wpg(ipg) * det(jac);
        end
    end
    eta_el(iele) = sqrt( e2_el(iele) / (e2_el(iele) + U2_el(iele)) );
    
end

etaG = sqrt( sum(e2_el) / (sum(e2_el) + sum(U2_el)) )


%% Configuracion deformada
D = (reshape(D,nDofNod,[]))';
nodePosition = nodes + D(:,1:2);


%Graficacion
limites = [min(min(stressExtra(:,:,2))),max(max(stressExtra(:,:,2)))];

figure(1)
bandplot(elements,nodePosition,stress(:,:,2),limites,'k');
title('Tensiones en los nodos.')

figure(2)
bandplot(elements,nodePosition,stressExtra(:,:,2),limites,'k');
title('Tensiones extrapoladas de los puntos de superconvergencia.')

figure(3)
scalarbandplot(elements,nodePosition,eta_el,[],'k',[],'flat');


