% ELEMENTOS QUADRIL?TEROS
clear
close all
format short g

eleType = 'Q4'; %'Q4', 'Q8', 'Q9'

% Discretizacion
nodes = load('nodQ4_5.txt');
nodes = nodes(:,2:3);


% Q4
%elements = [1 2 3 4];   % Matriz de conectividades - ajustar de acuerdo a n?mero de nodos por elemento 
elements = load('eleQ4_5.txt');
elements = elements(:,2:(end-5));
 
% Q8
%  elements = load('eleQ4_2.txt');
%  elements = elements(:,2:(end-1));% Matriz de conectividades - ajustar de acuerdo a n?mero de nodos por elemento
% Q9
% elements = [1 2 3 4 5 6 7 8 9];   % Matriz de conectividades - ajustar de acuerdo a n?mero de nodos por elemento


nDofNod = 2;                    % N?mero de grados de libertad por nodo
nel = size(elements,1);         % N?mero de elementos
nNod = size(nodes,1);           % N?mero de nodos
nNodEle = size(elements,2);     % N?mero de nodos por elemento
nDofTot = nDofNod*nNod;         % N?mero de grados de libertad
nDims = size(nodes,2);          % N?mero de dimensiones del problema

% Propiedades del Material
t = 1; % Espesor
L = 100;
c = 10;

E = 1000;
NU = 0.25;
C = E/(1 - NU^2)*[ 1.0     NU         0.0
                    NU    1.0         0.0
                   0.0    0.0     (1 - NU)/2 ];
               
G = E/(2*(1+NU));
I = (2*c)^3*t/12;
               
R = zeros(nNod,nDofNod);
P = 80;
sigmaX = @(x,y) -(3/2)*P*x.*y/10^3;
sigmaY = @(x,y) 0;
tauXY = @(x,y) -(3/4)*(P/c)*(1-(y/c).^2);

uExacto = @(x,y) -(P*(x.^2 - L^2).*y)/(2*E*I) - NU*P.*y.*(y.^2-c^2)/(6*E*I) + P.*y.*(y.^2-c^2)/(6*G*I);
vExacto = @(x,y) (NU*P.*x.*(y.^2)/(2*E*I)) + (P*(x.^3-L^3)/(6*E*I))-((P*L^2/(2*E*I))+(NU*P*c^2/(6*E*I))+(P*c^2/(3*G*I)))*(x-L);

syms x y

uFun = -(P*(x.^2 - L^2)*y)/(2*E*I) - NU*P*y*(y^2-c^2)/(6*E*I) + P*y*(y^2-c^2)/(6*G*I);
vFun = (NU*P*x*(y^2)/(2*E*I)) + (P*(x^3-L^3)/(6*E*I))-((P*L^2/(2*E*I))+(NU*P*c^2/(6*E*I))+(P*c^2/(3*G*I)))*(x-L);
epsX = diff(uFun,x);
epsY = diff(vFun,y);
gama = diff(uFun,y) + diff(vFun,x);

epsX = matlabFunction(epsX);
epsY = matlabFunction(epsY);
gama = matlabFunction(gama);

rsInt = ones(1,2)*2; %*1,*2,*3
[wpg, upg, npg] = gauss(rsInt);

% Matriz de rigidez
K = zeros(nDofTot);
for iele = 1:nel
    Ke = zeros(nDofNod*nNodEle);
    nodesEle = nodes(elements(iele,:),:);
    for ipg = 1:npg
        % Punto de Gauss
        ksi = upg(ipg,1);
        eta = upg(ipg,2);  
        % Derivadas de las funciones de forma respecto de ksi, eta
        dN = shapefunsder([ksi eta],eleType);  
        % Derivadas de x,y, respecto de ksi, eta
        jac = dN*nodesEle;                      
        % Derivadas de las funciones de forma respecto de x,y.
        dNxy = jac\dN;          % dNxy = inv(jac)*dN
        
        B = zeros(size(C,2),nDofNod*nNodEle);
        B(1,1:2:2*nNodEle-1) = dNxy(1,:);
        B(2,2:2:2*nNodEle) = dNxy(2,:);
        B(3,1:2:2*nNodEle-1) = dNxy(2,:);
        B(3,2:2:2*nNodEle) = dNxy(1,:); 

        Ke = Ke + t*B'*C*B*wpg(ipg)*det(jac);
    end
    eleDofs = node2dof(elements(iele,:),nDofNod);
    K(eleDofs,eleDofs) = K(eleDofs,eleDofs) + Ke;  
end



for iele = 1:nel
    Ry = zeros(2,1);
    Rx = zeros(2,1);
    nodesEle = nodes(elements(iele,:),:);
    if(any(nodesEle(:,1) == 0))
        X = nodesEle(nodesEle(:,1) == 0,1);
        Y = nodesEle(nodesEle(:,1) == 0,2);
    for ipg = 1:2
        % Punto de Gauss
        ksi = -1;
        eta = upg(ipg,1);  
        % Derivadas de las funciones de forma respecto de ksi, eta
        dN = shapefunsder([ksi eta],eleType);  
        % Derivadas de x,y, respecto de ksi, eta
        jac = dN*nodesEle;                      
        % Derivadas de las funciones de forma respecto de x,y.
        dNxy = jac\dN;          % dNxy = inv(jac)*dN
        N = shapefuns([ksi eta],eleType);
        N = N(nodesEle(:,1) == 0);
  
        Ry = Ry + N'*N*t*jac(2,2)*(-tauXY(X,Y));
    end
    eleDofs = node2dof(elements(iele,:),nDofNod);
    R(elements(iele,nodesEle(:,1) == 0),2) = R(elements(iele,nodesEle(:,1) == 0),2) +Ry;
    end
    
     if(any(nodesEle(:,1) == L))
        X = nodesEle(nodesEle(:,1) == L,1);
        Y = nodesEle(nodesEle(:,1) == L,2);
    for ipg = 1:2
        % Punto de Gauss
        ksi = 1;
        eta = upg(ipg,1);  
        % Derivadas de las funciones de forma respecto de ksi, eta
        dN = shapefunsder([ksi eta],eleType);  
        % Derivadas de x,y, respecto de ksi, eta
        jac = dN*nodesEle;                      
        % Derivadas de las funciones de forma respecto de x,y.
        dNxy = jac\dN;          % dNxy = inv(jac)*dN
        N = shapefuns([ksi eta],eleType);
        N = N(nodesEle(:,1) == L);
        Rx = Rx + N'*N*t*jac(2,2)*(sigmaX(X,Y));  %% Ac? no va el menos. Si lo pongo da cualquier cosa
        Ry = Ry + N'*N*t*jac(2,2)*tauXY(X,Y);
    end
    
    R(elements(iele,nodesEle(:,1) == L),1) = R(elements(iele,nodesEle(:,1) == L),1) + Rx;
    R(elements(iele,nodesEle(:,1) == L),2) = R(elements(iele,nodesEle(:,1) == L),2) + Ry;
    end
end


%% Condiciones de Borde

bc = zeros(nNod, nDofNod);
bc(nodes(:,1) == L & (nodes(:,2) == c | nodes(:,2) == -c | nodes(:,2) == 0) ,1) = true;
bc(nodes(:,1) == L & nodes(:,2) == 0,2) = true;

free = ~bc;
free = reshape(free',[],1);
Rr = reshape(R',[],1);
Rr = Rr(free);

Dr = K(free,free)\Rr;

D = zeros(nDofTot,1);
D(free) = Dr;

stress = zeros(nel,nNodEle,3);
strain = zeros(nel,nNodEle,3);

switch eleType
    case 'Q4'
        uNod = [ -1 -1
                  1 -1
                  1  1
                 -1  1 ];
    case 'Q8'
        uNod = [ -1 -1
                  0 -1
                  1 -1
                  1  0
                  1  1
                  0  1
                 -1  1 
                 -1  0];
end
for iele = 1:nel
    nodesEle = nodes(elements(iele,:),:);
    for inode = 1:nNodEle
        % Punto de Gauss
        ksi = uNod(inode,1);
        eta = uNod(inode,2);  
        % Derivadas de las funciones de forma respecto de ksi, eta
         dN = shapefunsder([ksi eta],eleType);  
        % Derivadas de x,y, respecto de ksi, eta
        jac = dN*nodesEle;                      
        % Derivadas de las funciones de forma respecto de x,y.
        dNxy = jac\dN;          % dNxy = inv(jac)*dN
        
        B = zeros(size(C,2),nDofNod*nNodEle);
        B(1,1:2:2*nNodEle-1) = dNxy(1,:);
        B(2,2:2:2*nNodEle) = dNxy(2,:);
        B(3,1:2:2*nNodEle-1) = dNxy(2,:);
        B(3,2:2:2*nNodEle) = dNxy(1,:);
        
        eleDofs = node2dof(elements(iele,:),nDofNod);
        eleDofs = reshape(eleDofs',[],1);
        
        stress(iele,inode,:) = C*B*D(eleDofs);
  
    end
end

for iele = 1:nel
    nodesEle = nodes(elements(iele,:),:);
    for inode = 1:nNodEle
        % Punto de Gauss
        ksi = uNod(inode,1);
        eta = uNod(inode,2);  
        % Derivadas de las funciones de forma respecto de ksi, eta
         dN = shapefunsder([ksi eta],eleType);  
        % Derivadas de x,y, respecto de ksi, eta
        jac = dN*nodesEle;                      
        % Derivadas de las funciones de forma respecto de x,y.
        dNxy = jac\dN;          % dNxy = inv(jac)*dN
        
        B = zeros(size(C,2),nDofNod*nNodEle);
        B(1,1:2:2*nNodEle-1) = dNxy(1,:);
        B(2,2:2:2*nNodEle) = dNxy(2,:);
        B(3,1:2:2*nNodEle-1) = dNxy(2,:);
        B(3,2:2:2*nNodEle) = dNxy(1,:);
        
        eleDofs = node2dof(elements(iele,:),nDofNod);
        eleDofs = reshape(eleDofs',[],1);
        
        strain(iele,inode,:) = B*D(eleDofs);
    end
end


eleError =  zeros(nel,1);
nodeErrorVector = zeros(nel,nNodEle);

for iele = 1:nel
    nodesEle = nodes(elements(iele,:),:);
    
    for inode = 1:nNodEle
        aproxStrain = zeros(3,1);
        aproxStrain(1) = strain(iele,inode,1);
        aproxStrain(2) = strain(iele,inode,2);
        aproxStrain(3) = strain(iele,inode,3);
      
        X = nodesEle(inode,1);
        Y = nodesEle(inode,2);
        exactStrain = [epsX(X,Y); epsY(X,Y); gama(X,Y)];  %% Aca van las funciones que no calcule todavia
        nodeError = 0;
        for ipg = 1:npg
        % Punto de Gauss
        ksi = upg(ipg,1);
        eta = upg(ipg,2);  
        % Derivadas de las funciones de forma respecto de ksi, eta
         dN = shapefunsder([ksi eta],eleType);  
        % Derivadas de x,y, respecto de ksi, eta
        jac = dN*nodesEle;                      
        % Derivadas de las funciones de forma respecto de x,y.
        dNxy = jac\dN;          % dNxy = inv(jac)*dN
   
        N = shapefuns([ksi eta],eleType);
        N = N(inode);
        
        nodeError = nodeError +(aproxStrain-exactStrain)'*C*(aproxStrain-exactStrain);
     
        end
        nodeErrorVector(iele,inode) = nodeError;
    end
    eleError(iele) = eleError(iele) + nodeError;
    
end

eleErrorPercent = 100*eleError;

bad = eleErrorPercent > 5;
badEle = 1:nel;
badEle = badEle(bad);





D = reshape(D,2,[])';
nodePosition = nodes + D;
exactSolution = nodes + [uExacto(nodes(:,1),nodes(:,2)) ,vExacto(nodes(:,1),nodes(:,2))];
figure
meshplot(elements,nodes,'b')
meshplot(elements,nodePosition,'r')

meshplot(elements,exactSolution,'g')
figure
bandplot(elements,nodePosition,stress(:,:,1),[],'k');
figure
bandplot(elements,nodePosition,100*nodeErrorVector(:,:),[],'k');



