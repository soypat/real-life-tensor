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
t=1;
E = 200e3;
NU = 0.3;
 C = E/((1 + NU)*(1 - 2*NU))*[ 1 - NU      NU            0.0
                                      NU       1 - NU         0.0
                                      0.0        0.0     (1 - 2*NU)/2 ];
                                  
%% Ensamble
%% Gauss           
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
        y = upg(ipg,2);  %% las func de forma son iguales a las calculadas antes 
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

autovals = eig(K);          % 3 dan 0. Movimientos Rigidos
       