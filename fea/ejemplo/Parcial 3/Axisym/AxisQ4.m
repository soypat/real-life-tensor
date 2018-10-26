
clear
clc
close all
format short g

load('Nodos2x6Q4.mat');
load('Elementos2x6Q4.mat');

nDofNod = 2;                    % Numero de grados de libertad por nodo
nNodEle = size(elements,2);     % Numero de nodos por elemento
nel = size(elements,1);         % Numero de elementos
nNod = size(nodes,1);           % Numero de nodos
nDofTot = nDofNod*nNod;         % Numero de grados de libertad

bc = false(nNod,nDofNod);       % Matriz de condiciones de borde
bc(abs(nodes(:,2)-0) <= 1e-4 ,2) = true;
bc(abs(nodes(:,2)-2) <= 1e-4 ,2) = true;


r1 = 4; %mm
r2 = 10; %mm
d = 2; %mm

%% Gauss           
a   = 1/sqrt(3);
% Ubicaciones puntos de Gauss
upg = [ -a  -a
         a  -a
         a   a
        -a   a ];
% Numero de puntos de Gauss
npg = size(upg,1);
wpg = ones(4,1);

%% Matriz Constitutiva
E = 1;
NU = 0.3;
f = NU/(1 - NU);
g = (1 - 2*NU)/(2*(1 - NU));
C = [ 1 f f 0
      f 1 f 0
      f f 1 0
      0 0 0 g ] * (1 - NU)*E/((1 + NU)*(1 - 2*NU));
               
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
        dN = shapefunsder([ksi eta],'Q4');
        N = shapefuns([ksi eta],'Q4');
        % Derivadas de x,y, respecto de ksi, eta
        jac = dN*nodesEle;                      
        % Derivadas de las funciones de forma respecto de x,y.
        dNxy = jac\dN;          % dNxy = inv(jac)*dN
        
        r = N*nodesEle(:,1);
        
        B = zeros(size(C,2),nDofNod*nNodEle);
        B(1,1:2:7) = dNxy(1,:);
        B(2,1:2:7) = N/r;
        B(3,2:2:8) = dNxy(2,:);
        B(4,1:2:7) = dNxy(2,:);
        B(4,2:2:8) = dNxy(1,:);  

        Ke = Ke + B'*C*B*wpg(ipg)*det(jac)*r;
    end
    eleDofs = nodeDofs(elements(iele,:),:);
    eleDofs = reshape(eleDofs',[],1);
    K(eleDofs,eleDofs) = K(eleDofs,eleDofs) + Ke;  
end

%% Carga
p = 100e-1; % [N/mm^2]
% Ubicaciones puntos de Gauss
upg = [-a a];
% Numero de puntos de Gauss
npg = size(upg,2);
wpg = [1 1];
R = zeros(nNod,nDofNod);        % Vector de cargas
nodeDofs = reshape(1:nDofTot,nDofNod,nNod)';
for iele = 1:nel
    nodesEle = nodes(elements(iele,:),:);
    if abs(nodesEle(1,1)-r1) <= 1e-4
        for ipg = 1:npg
            % Punto de Gauss
            ksi = -1;
            eta = upg(ipg);  
            % Derivadas de las funciones de forma respecto de ksi, eta
            dN = shapefunsder([ksi eta],'Q4'); 
            % Derivadas de x,y, respecto de ksi, eta
            jac = dN*nodesEle;                     
            % Derivadas de las funciones de forma respecto de x,y.
            dNxy = jac\dN;          % dNxy = inv(jac)*dN
            N = shapefuns([ksi eta],'Q4');
            r = N*nodesEle(:,1);
            N = N([1 4]);
            R(elements(iele,[1 4]),1) = R(elements(iele,[1 4]),1) + ...
                N'*N*p*ones(2,1)*jac(2,2)*wpg(ipg)*r;
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

% D = (reshape(D,nDofNod,[]))';
% nodePosition = nodes + D(:,1:2);

% figure
% meshplot(elements,nodes,'b')
% hold on
% meshplot(elements,nodePosition,'r')

%% Recuperacion de tensiones en los nodos
stress = zeros(nel,nNodEle,4);
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
        N = shapefuns([ksi eta],'Q4'); 
        % Derivadas de x,y, respecto de ksi, eta
        jac = dN*nodesEle;                      
        % Derivadas de las funciones de forma respecto de x,y.
        dNxy = jac\dN;          % dNxy = inv(jac)*dN
        
        r = N*nodesEle(:,1);
        
        B = zeros(size(C,2),nDofNod*nNodEle);
        B(1,1:2:7) = dNxy(1,:);
        B(2,1:2:7) = N/r;
        B(3,2:2:8) = dNxy(2,:);
        B(4,1:2:7) = dNxy(2,:);
        B(4,2:2:8) = dNxy(1,:);  
        
        eleDofs = nodeDofs(elements(iele,:),:);
        eleDofs = reshape(eleDofs',[],1);
        stress(iele,inode,:) = C*B*D(eleDofs);
    end
end

%% tensiones en los puntos de superconvergencia
stressSuper = zeros(nel,nNodEle,4);
for iele = 1:nel
    nodesEle = nodes(elements(iele,:),:);
    ksi = 0;
    eta = 0; 
    dN = shapefunsder([ksi eta],'Q4');
    N = shapefuns([ksi eta],'Q4');
    for inode = 1:nNodEle
        jac = dN*nodesEle;                      
        % Derivadas de las funciones de forma respecto de x,y.
        dNxy = jac\dN;          % dNxy = inv(jac)*dN
        
        r = N*nodesEle(:,1);
        
        B = zeros(size(C,2),nDofNod*nNodEle);
        B(1,1:2:7) = dNxy(1,:);
        B(2,1:2:7) = N/r;
        B(3,2:2:8) = dNxy(2,:);
        B(4,1:2:7) = dNxy(2,:);
        B(4,2:2:8) = dNxy(1,:); 
        
        eleDofs = nodeDofs(elements(iele,:),:);
        eleDofs = reshape(eleDofs',[],1);
        stressSuper(iele,inode,:) = C*B*D(eleDofs);
    end
end

%% promediado de tensiones en los nodos
avgStress = zeros(nNod,4);
for inode = 1:nNod
    [I,J] = find(elements == inode);
    nShare = length(I);
    for ishare = 1:nShare
        avgStress(inode,:) = avgStress(inode,:) + squeeze(stressSuper(I(ishare),J(ishare),:))';
    end
    avgStress(inode,:) = avgStress(inode,:) / nShare;
end

%% Configuracion deformada
D = (reshape(D,nDofNod,[]))';
nodePosition = nodes + D(:,1:2);

%Graficacion
% limites = [min(min(stressExtra(:,:,2))),max(max(stressExtra(:,:,2)))];

figure(1)
bandplot(elements,nodePosition,stress(:,:,1),[],'k');
title('Tensiones en el centro de cada elemento.')

figure(2)
scalarbandplot(elements,nodePosition,avgStress(:,1),[],'k',[],'interp');
title('Tensiones promediadas.')

% figure(3)
% scalarbandplot(elements,nodePosition,eta_el,[],'k',[],'flat');

%% Error
% %% Gauss           
% [wpg, upg, npg] = gauss([2 2]);
% 
% %% calculo de eta_el, e2, U2
% invC = C\eye(3);
% eta_el = zeros(nel,1);
% e2_el = zeros(nel,1);
% U2_el = zeros(nel,1);
% for iele = 1:nel
%     nodesEle = nodes(elements(iele,:),:);
%     for inode = 1:nNodEle
%         for ipg = 1:npg
%             % Punto de Gauss
%             ksi = upg(ipg,1);
%             eta = upg(ipg,2);
% 
%             % Derivadas de las funciones de forma respecto de ksi, eta
%             dN = shapefunsder([ksi eta],'Q8');
%             % Derivadas de x,y, respecto de ksi, eta
%             jac = dN*nodesEle;                      
%             % Derivadas de las funciones de forma respecto de x,y.
%             dNxy = jac\dN;          % dNxy = inv(jac)*dN
% 
%             % funciones de forma
%             N = shapefuns([ksi eta],'Q8');
%             eleStress =  squeeze(stress(iele,inode,:));            % tensiones "directas"
%             starStress = squeeze(stressExtra(iele,inode,:)); % tensiones mejoradas
% 
%             e2_el(iele) = e2_el(iele) + (starStress - eleStress)' * ... 
%                         invC * (starStress - eleStress) * wpg(ipg) * det(jac);
% 
%             U2_el(iele) = U2_el(iele) + eleStress' * invC * eleStress * ...
%                           wpg(ipg) * det(jac);
%         end
%     end
%     eta_el(iele) = sqrt( e2_el(iele) / (e2_el(iele) + U2_el(iele)) );
%     
% end
% 
% etaG = sqrt( sum(e2_el) / (sum(e2_el) + sum(U2_el)) )
