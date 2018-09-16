clear 
clc
close all


%% Discretización
N = 100;
condNbr = zeros(N,2);
wt = 10;
Ro = linspace(wt,10000*wt,N);
H = 1;

for i = 1:N
R1 = Ro(i) + wt;

nodes = [ Ro(i) 0
          R1    0
          R1    H
          Ro(i) H ];

elements = [1  2  3  4];      %Matriz de conectividades


eleType = 'Q4';
nDofNod = 2;                    % Número de grados de libertad por nodo
switch eleType
    case 'Q4'
        nNodEle = 4;            % Número de nodos por elemento
end
nel = size(elements,1);         % Número de elementos
nNod = size(nodes,1);           % Número de nodos
nDofTot = nDofNod*nNod;         % Número de grados de libertad

bc = false(nNod,nDofNod);       % Matriz de condiciones de borde
bc(1,2) = true;


%% Propiedades del material
E = 1;
NU = 0.33;
%%
meshplot(elements,nodes,'b')

%% Matriz Constitutiva

f = NU/(1 - NU);
g = (1 - 2*NU)/(2*(1 - NU));
C = [ 1 f f 0
      f 1 f 0
      f f 1 0
      0 0 0 g ] * (1 - NU)*E/((1 + NU)*(1 - 2*NU));

%% Matriz de rigidez
% Gauss           
rsInt = 2*ones(1,2);
[wpg, upg, npg] = gauss(rsInt);
% Funciones de forma y derivadas
Ni = shapefuns   (upg,eleType);
dN = shapefunsder(upg,eleType);

K = zeros(nDofTot);
for iele = 1:nel
    Ke = zeros(nDofNod*nNodEle);
    nodesEle = nodes(elements(iele,:),:);
    for ipg = 1:npg
        % Derivadas de x,y, respecto de ksi, eta
        jac = dN(:,:,ipg)*nodesEle;                      
        % Derivadas de las funciones de forma respecto de x,y.
        dNxy = jac\dN(:,:,ipg);          % dNxy = inv(jac)*dN(:,:,ipg)

        r = Ni(:,:,ipg)*nodesEle(:,1);
        
        B = zeros(size(C,2),nDofNod*nNodEle);
        B(1,1:2:7) = dNxy(1,:);
        B(2,1:2:7) = Ni(:,:,ipg)/r;
        B(3,2:2:8) = dNxy(2,:);
        B(4,1:2:7) = dNxy(2,:);
        B(4,2:2:8) = dNxy(1,:);        
        
        Ke = Ke + B'*C*B*r*det(jac)*wpg(ipg);
    end
    eleDofs = node2dof(elements(iele,:),nDofNod);
    K(eleDofs,eleDofs) = K(eleDofs,eleDofs) + Ke;  
end
%% Reduccion Matriz
isFixed = reshape(bc',[],1);
isFree = ~isFixed;

condNbr(i,2) = cond(K(isFree,isFree));
condNbr(i,1) = ( (Ro(i) + R1)/2 ) / wt;
end

plot(condNbr(:,1),condNbr(:,2),'-*')
grid