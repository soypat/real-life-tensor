
clear
clc
close all
format short g

load('Nodos2x4Q8.mat');
load('Elementos2x4Q8.mat');

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
[wpg, upg, npg] = gauss([3 3]);

%% Matriz Constitutiva (plane stress)
E = 1000;
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
        dN = shapefunsder([ksi eta],'Q8');
        N = shapefuns([ksi eta],'Q8');
        % Derivadas de x,y, respecto de ksi, eta
        jac = dN*nodesEle;                      
        % Derivadas de las funciones de forma respecto de x,y.
        dNxy = jac\dN;          % dNxy = inv(jac)*dN
        
        r = N*nodesEle(:,1);
        
        B = zeros(size(C,2),nDofNod*nNodEle);
        B(1,1:2:nNodEle*nDofNod-1) = dNxy(1,:);
        B(2,1:2:nNodEle*nDofNod-1) = N/r;
        B(3,2:2:nNodEle*nDofNod)   = dNxy(2,:);
        B(4,1:2:nNodEle*nDofNod-1) = dNxy(2,:);
        B(4,2:2:nNodEle*nDofNod)   = dNxy(1,:);  

        Ke = Ke + B'*C*B*wpg(ipg)*det(jac)*r;
    end
    eleDofs = nodeDofs(elements(iele,:),:);
    eleDofs = reshape(eleDofs',[],1);
    K(eleDofs,eleDofs) = K(eleDofs,eleDofs) + Ke;  
end

%% Carga
omega = 0.5; % [rad/s]
rho = 3; % [kg/m^3]
o2r = omega^2*rho;
R = zeros(nNod,nDofNod);        % Vector de cargas
[wpg, upg, npg] = gauss([2 2]);
% upg = upg([1 3 4 2]);
for iele = 1:nel
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
        N = shapefuns([ksi eta],'Q8');
        r = N*nodesEle(:,1);
        R(elements(iele,:),1) = R(elements(iele,:),1) + ...
            N'*o2r*r^2*det(jac)*wpg(ipg);
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
% 
% figure
% meshplot(elements,nodes,'b')
% hold on
% meshplot(elements,nodePosition,'r')

%% Recuperacion de tensiones en los nodos
stress = zeros(nel,nNodEle,4);
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
        N = shapefuns([ksi eta],'Q8');
        % Derivadas de x,y, respecto de ksi, eta
        jac = dN*nodesEle;                      
        % Derivadas de las funciones de forma respecto de x,y.
        dNxy = jac\dN;          % dNxy = inv(jac)*dN
        
        r = N*nodesEle(:,1);
        
        B = zeros(size(C,2),nDofNod*nNodEle);
        B(1,1:2:nNodEle*nDofNod-1) = dNxy(1,:);
        B(2,1:2:nNodEle*nDofNod-1) = N/r;
        B(3,2:2:nNodEle*nDofNod)   = dNxy(2,:);
        B(4,1:2:nNodEle*nDofNod-1) = dNxy(2,:);
        B(4,2:2:nNodEle*nDofNod)   = dNxy(1,:); 
        
        eleDofs = nodeDofs(elements(iele,:),:);
        eleDofs = reshape(eleDofs',[],1);
        stress(iele,inode,:) = C*B*D(eleDofs);
    end
end

%% tensiones en los puntos de superconvergencia
stressSuper = zeros(nel,4,4);
uNod = [-1 -1
         1 -1
         1  1
        -1  1]/sqrt(3);
for iele = 1:nel
    nodesEle = nodes(elements(iele,:),:); 
    for inode = 1:4
        ksi = uNod(inode,1);
        eta = uNod(inode,2);
        dN = shapefunsder([ksi eta],'Q8');
        N = shapefuns([ksi eta],'Q8');
        % Derivadas de x,y, respecto de ksi, eta
        jac = dN*nodesEle;                      
        % Derivadas de las funciones de forma respecto de x,y.
        dNxy = jac\dN;          % dNxy = inv(jac)*dN

        r = N*nodesEle(:,1);

        B = zeros(size(C,2),nDofNod*nNodEle);
        B(1,1:2:nNodEle*nDofNod-1) = dNxy(1,:);
        B(2,1:2:nNodEle*nDofNod-1) = N/r;
        B(3,2:2:nNodEle*nDofNod)   = dNxy(2,:);
        B(4,1:2:nNodEle*nDofNod-1) = dNxy(2,:);
        B(4,2:2:nNodEle*nDofNod)   = dNxy(1,:);
        
        eleDofs = nodeDofs(elements(iele,:),:);
        eleDofs = reshape(eleDofs',[],1);
        stressSuper(iele,inode,:) = C*B*D(eleDofs);
    end
end

% Extrapolación a los Nodos
rsExt = [-1  -1
          1  -1
          1   1
         -1   1
          0  -1
          1   0
          0   1
         -1   0]*sqrt(3);
stressExtra = zeros(nel,nNodEle,4);
for iele = 1:nel
    for inode = 1:nNodEle
        rr = rsExt(inode,1);
        s = rsExt(inode,2);
        
        N = shapefuns([rr s],'Q4');
        
        stressExtra(iele,inode,:) = N * squeeze(stressSuper(iele,:,:));
    end
end


%% promediado de tensiones en los nodos
avgStress = zeros(nNod,4);
for inode = 1:nNod
    [I,J] = find(elements == inode);
    nShare = length(I);
    for ishare = 1:nShare
        avgStress(inode,:) = avgStress(inode,:) + squeeze(stressExtra(I(ishare),J(ishare),:))';
    end
    avgStress(inode,:) = avgStress(inode,:) / nShare;
end

%% Configuracion deformada
D = (reshape(D,nDofNod,[]))';
nodePosition = nodes + D(:,1:2);

%Graficacion
% limites = [min(min(stressExtra(:,:,2))),max(max(stressExtra(:,:,2)))];

figure(1)
bandplot(elements,nodePosition,stress(:,:,2),[],'k');
title('Tensiones en los nodos.')

figure(2)
bandplot(elements,nodePosition,stressExtra(:,:,2),[],'k');
title('Tensiones extrapoladas.')

figure(3)
scalarbandplot(elements,nodePosition,avgStress(:,2),[],'k',[],'interp');
title('Tensiones promediadas en los nodos.')

% figure(3)
% scalarbandplot(elements,nodePosition,eta_el,[],'k',[],'flat');

vm = sqrt(((stressExtra(:,:,1)-stressExtra(:,:,2)).^2 + ...
    (stressExtra(:,:,2)-stressExtra(:,:,3)).^2 + ...
    (stressExtra(:,:,3)-stressExtra(:,:,1)).^2 + ...
    6*stressExtra(:,:,4).^2)/2);

figure
bandplot(elements,nodePosition,vm,[],'k');
title('VM')

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
