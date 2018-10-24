%% An application with high stress gradient
% close all; clear; clc
%% Datos
% 1.Acero 2.Aluminio
graph=true;
E = [200 70]*1e3; % [MPa]
NU = [.29 .33];
alfa = [12 24]*1e-6; %[1/C]
eleT = 'Q8';
dT = 1000; % [C]
t = 20; % [mm]
%% Nodos y elementos

elementos = load('elem_bimetallic.txt');
elementos = elementos(:,2:9);

nodes = load('nod_bimetallic.txt'); % [mm]
nodes = nodes(:,2:3);
C1 = E(1)/(1 - NU(1)^2)*[ 1.0     NU(1)       0.0
                                NU(1)     1.0        0.0
                                 0.0      0.0    (1 - NU(1))/2 ];
C2 = E(2)/(1 - NU(2)^2)*[ 1.0     NU(2)       0.0
                                NU(2)     1.0        0.0
                                 0.0      0.0    (1 - NU(2))/2 ];
                             
                             
Ndofpornodo = 2;                    % Número de grados de libertad por nodo
Nnodporelem = 8;                    % Número de nodos por elemento
Nelem = size(elementos,1);         % Número de elementos
Nnod = size(nodes,1);           % Número de nodos
doftot = Ndofpornodo*Nnod;         % Número de grados de libertad
DOF=reshape(1:doftot,2,[])';
bc = false(Nnod,Ndofpornodo);       % Matriz de condiciones de borde
R = zeros(Nnod,Ndofpornodo);        % Vector de cargas
if graph==true
figure(1)
Meshplot(elementos,nodes,bc,'k',1)
end
%% Gauss para regla de 3x3 (numeración de nodos como figura 6.3-3 pág 212)          
a   = sqrt(0.6);
% Ubicaciones puntos de Gauss
upg = [ -a  -a
        -a   0
        -a   a
         0  -a    
         0   0
         0   a 
         a  -a
         a   0
         a   a ];
% Número de puntos de Gauss
npg = size(upg,1);
% wpg=5/9*ones(npg,1);wpg(5)=8/9; FORMA INCORRECTA
W = [5/9 8/9 5/9];
wpg = reshape(W'*W,1,[]);

%% Matriz de rigidez y deformaciones térmicas
K = zeros(doftot);
nodeDofs = reshape(1:doftot,Ndofpornodo,Nnod)';
epsilon0 = dT*alfa;
P = zeros(doftot,1);
for iele = 1:Nelem
    index=elementos(iele,:);
    meindof = reshape(nodeDofs(index,:)',[],1);
    Ke = zeros(Ndofpornodo*Nnodporelem);
    Pe = zeros(Ndofpornodo*Nnodporelem,1);
    nodesEle = nodes(index,:);
    if nodesEle(4,2) == 60
        C = C1;
        eps = epsilon0(1);
        epsilon = [eps;eps;0];
    else
        C = C2;
        eps = epsilon0(2);
        epsilon = [eps;eps;0];
    end
    for ipg = 1:npg
        ksi = upg(ipg,1);
        eta = upg(ipg,2);  
        dN = shapefunsder([ksi eta],eleT); 
        jac = dN*nodesEle;                      
        dNxy = jac\dN;
        
        B = zeros(size(C,2),Ndofpornodo*Nnodporelem);
        B(1,1:2:15) = dNxy(1,:);
        B(2,2:2:16) = dNxy(2,:);
        B(3,1:2:15) = dNxy(2,:);
        B(3,2:2:16) = dNxy(1,:); 

        Ke = Ke + B'*C*B*wpg(ipg)*det(jac)*t;
        Pe = Pe + B'*C*epsilon*wpg(ipg)*det(jac)*t;
    end
    K(meindof,meindof) = K(meindof,meindof) + Ke;
    P(meindof) = P(meindof) + Pe;
end


%% BC
bc(1,:) = true;
bc(2,2) = true;

%% Reduccion Matriz
isFixed = reshape(bc',[],1);
isFree = ~isFixed;

Rr = reshape(P,[],1);

% Solver
Dr = K(isFree,isFree)\Rr(isFree);

% Reconstrución
D = zeros(doftot,1);
D(isFree) = D(isFree) + Dr;

%% Tensiones en los nodos
stress = zeros(Nelem,Nnodporelem,3);
unod = [ -1 -1
          1 -1
          1  1
         -1  1
          0 -1
          1  0
          0  1
         -1  0];
for iele = 1:Nelem
    index=elementos(iele,:);
    nodesEle = nodes(index,:);
    if nodesEle(4,2) == 60
        C = C1;
        eps = epsilon0(1);
        epsilon = [eps;eps;0];
    else
        C = C2;
        eps = epsilon0(2);
        epsilon = [eps;eps;0];
    end
    for in = 1:Nnodporelem
        % Punto de Gauss
        ksi = unod(in,1);
        eta = unod(in,2);  
        % Derivadas de las funciones de forma respecto de ksi, eta
        dNke = shapefunsder([ksi eta],eleT);
        % Derivadas de x,y, respecto de ksi, eta
        J = dNke*nodesEle;                      
        % Derivadas de las funciones de forma respecto de x,y.
        dNxy = J\dNke;
        B = zeros(size(C,2),Ndofpornodo*Nnodporelem);
        B(1,1:2:15) = dNxy(1,:);
        B(2,2:2:16) = dNxy(2,:);
        B(3,1:2:15) = dNxy(2,:);
        B(3,2:2:16) = dNxy(1,:);
        dofs = reshape(nodeDofs(index,:)',[],1);
        stress(iele,in,:) = C*(B*D(dofs)-epsilon);
    end
end

%% Tensiones en los puntos de superconvergencia
RSextrapolation = unod*sqrt(3);
a   = 1/sqrt(3);
% Ubicaciones puntos de Gauss
upg = [ -a  -a
        -a   a
         a  -a
         a   a ];    
% Número de puntos de Gauss
npg = size(upg,1);
wpg = ones(npg,1);
stress2 = zeros(Nelem,Nnodporelem,3);
for iele = 1:Nelem
    index=elementos(iele,:);
    nodesEle = nodes(index,:);
    stressgauss = zeros(npg,3);
    if nodesEle(4,2) == 60
        C = C1;
        eps = epsilon0(1);
        epsilon = [eps;eps;0];
    else
        C = C2;
        eps = epsilon0(2);
        epsilon = [eps;eps;0];
    end
    for ipg = 1:npg
        % Punto de Gauss
        ksi = upg(ipg,1);
        eta = upg(ipg,2);  
        % Derivadas de las funciones de forma respecto de ksi, eta
        dNke = shapefunsder([ksi eta],eleT);
        % Derivadas de x,y, respecto de ksi, eta
        J = dNke*nodesEle;                      
        % Derivadas de las funciones de forma respecto de x,y.
        dNxy = J\dNke;
        B = zeros(size(C,2),Ndofpornodo*Nnodporelem);
        B(1,1:2:15) = dNxy(1,:);
        B(2,2:2:16) = dNxy(2,:);
        B(3,1:2:15) = dNxy(2,:);
        B(3,2:2:16) = dNxy(1,:);
        dofs = reshape(nodeDofs(index,:)',[],1);
        stressgauss(ipg,:) = (C*(B*D(dofs)-epsilon))';
    end
    StressElem = zeros(Nnodporelem,3);
    for iNod = 1:8
        r = RSextrapolation(iNod,1);
        s = RSextrapolation(iNod,2);
        N = shapefuns([r s],'Q4');
        StressElem(iNod,:) = N*stressgauss;
    end
    stress2(iele,:,:) = StressElem;
end

%% Configuracion deformada
D = (reshape(D,Ndofpornodo,[]))';
nodePosition = nodes + D(:,1:2)*10;

%Gráfico
if graph==true
    figure(2)
    bandplot(elementos,nodePosition,stress2(:,:,1),[],'k');
end

































