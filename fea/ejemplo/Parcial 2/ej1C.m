clc
clear

nodes = [ -1 1
           1 1
           0 0
          -1 -1
           1 -1]*1e3;
       
elements = [ 4 5 2 1 3];

nDofNod = 2;
nNodEle = 5;
nNod = 5;
nEle = 1;
nDofTot = nNod*nDofNod;

%% Matriz Constitutiva (plane strain)

% Trabajo en N, mm, MPa
t=1;
E = 200e3;
NU = 0.3;
C = E/((1 + NU)*(1 - 2*NU))*[ 1 - NU      NU            0.0
                                      NU       1 - NU         0.0
                                      0.0        0.0     (1 - 2*NU)/2 ];
                                  
%% Ensamble
%% Gauss           2x2
a   = 1/sqrt(3);
% Ubicaciones puntos de Gauss
upg = [ -a  -a
         a  -a
         a   a
        -a   a ];    
% Número de puntos de Gauss
npg = size(upg,1);
wpg = ones(npg,1);

%% Matriz de rigidez (Para elementos Q4 con regla de Gauss de 2x2)
K = zeros(nDofTot);
nodeDofs = reshape(1:nDofTot,nDofNod,nNod)';

for iele = 1:nEle
    Ke = zeros(nDofNod*nNodEle);
    nodesEle = nodes(elements(iele,:),:);
    for ipg = 1:npg
        % Punto de Gauss
        x = upg(ipg,1);
        y = upg(ipg,2);  %% Es indentico al cuadrado de isoparametricas en x e y
        % Derivadas de las funciones de forma respecto de ksi, eta
        dNx = [ x/2 + y/4 - 1/4, x/2 - y/4 + 1/4, x/2 + y/4 + 1/4, x/2 - y/4 - 1/4, -2*x];
        dNy = [ x/4 - 1/4, - x/4 - 1/4, x/4 + 1/4, 1/4 - x/4, 0];
        dN = [dNx
              dNy];
        % Derivadas de x,y, respecto de ksi, eta
        jac = dN*nodesEle;                      
        % Derivadas de las funciones de forma respecto de x,y.
        dNxy = jac\dN;          % dNxy = inv(jac)*dN
        
        B = zeros(size(C,2),nDofNod*nNodEle);
        B(1,1:2:9) = dNxy(1,:);
        B(2,2:2:10) = dNxy(2,:);
        B(3,1:2:9) = dNxy(2,:);
        B(3,2:2:10) = dNxy(1,:); 

        Ke = Ke + B'*C*B*wpg(ipg)*det(jac);
    end
    eleDofs = nodeDofs(elements(iele,:),:);
    eleDofs = reshape(eleDofs',[],1);
    K(eleDofs,eleDofs) = K(eleDofs,eleDofs) + Ke;  
end

%% Condiciones de borde y cargas

bc = false(nNod,nDofNod);       % Matriz de condiciones de borde

R = zeros(nNod,nDofNod);        % Vector de cargas

bc(4,:) = true;
bc(5,2) = true;
%N =[ (x*y)/4 - y/4 - x/4 + x^2/4, x/4 - y/4 - (x*y)/4 + x^2/4, x/4 + y/4 + (x*y)/4 + x^2/4, y/4 - x/4 - (x*y)/4 + x^2/4, 1 - x^2]
syms x real
q = 2e-3;
qy =  q*x;
y = 1e3;
N4 =  y/4 - x/4 - (x*y)/4 + x^2/4 ;
N3 =  x/4 + y/4 + (x*y)/4 + x^2/4 ;
int1 =  qy*N4;
int2 =  qy*N3;
R(1,2) = int(int1,-1,1);
R(2,2) = int(int2,-1,1);

%% Reduccion Matriz
isFixed = reshape(bc',[],1);
isFree = ~isFixed;

Rr = reshape(R',[],1);

% Solver
Dr = K(isFree,isFree)\Rr(isFree);

% Reconstrución
D = zeros(nDofTot,1);
D(isFree) = D(isFree) + Dr;


%% Recuperacion de tensiones en los 4 ptos de gauss

stress = zeros(npg,3,nEle);
vm = zeros(npg,1);
for iele = 1:nEle
    nodesEle = nodes(elements(iele,:),:);
    for ipg = 1:npg
        % Punto de Gauss
        x = upg(ipg,1);
        y = upg(ipg,2);  %% Es indentico al cuadrado de isoparametricas en x e y
        % Derivadas de las funciones de forma respecto de ksi, eta
        dNx = [ x/2 + y/4 - 1/4, x/2 - y/4 + 1/4, x/2 + y/4 + 1/4, x/2 - y/4 - 1/4, -2*x];
        dNy = [ x/4 - 1/4, - x/4 - 1/4, x/4 + 1/4, 1/4 - x/4, 0];
        dN = [dNx
              dNy];
        % Derivadas de x,y, respecto de ksi, eta
        jac = dN*nodesEle;                      
        % Derivadas de las funciones de forma respecto de x,y.
        dNxy = jac\dN;          % dNxy = inv(jac)*dN
        
        B = zeros(size(C,2),nDofNod*nNodEle);
        B(1,1:2:9) = dNxy(1,:);
        B(2,2:2:10) = dNxy(2,:);
        B(3,1:2:9) = dNxy(2,:);
        B(3,2:2:10) = dNxy(1,:); 
        eleDofs = nodeDofs(elements(iele,:),:);
        eleDofs = reshape(eleDofs',[],1);
        stress(ipg,:,iele) = C*B*D(eleDofs);
        sigmax = stress(ipg,1,iele);
        sigmay = stress(ipg,2,iele);
        tau = stress(ipg,3,iele);
        
        vm(ipg) = sqrt(sigmax^2 -sigmax*sigmay+sigmay^2 +3*tau^2);
      
    end
end


maximo = max(vm);

    
   







