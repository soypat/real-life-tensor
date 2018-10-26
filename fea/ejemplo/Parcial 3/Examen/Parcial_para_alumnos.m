%% Parcial 3

% Alumno: 

clear
clc
close all

%% Datos
t=1;%mm
E=1; %MPa
nu=1/5;

S = 1; % escala

eleType = 'Q4';

%Plane Strain
C = (E/((1+nu)*(1-2*nu)))*[1-nu    nu      0;
                             nu  1-nu      0;
                              0    0  0.5-nu];

%% Discretización

if strcmp(eleType,'Q4')
    nNodEle = 4;
    Nodos = [  0    0 %Matriz de coordenadas Q4
        50    0
        100    0
        0   50
        50   50
        100   50
        0  100
        50  100
        60  100
        80   80
        100   60];
    
    Elementos = [1  2  5  4;      %Matriz de conectividades Q4
        2  3  6  5;
        4  5  8  7;
        5 10  9  8;
        5  6 11 10];
        
elseif strcmp(eleType,'Q5')
    nNodEle = 5;
    
    Nodos = [  0    0 %Matriz de coordenadas Q5
        50    0
        100    0
        0   50
        50   50
        100   50
        0  100
        50  100
        60  100
        80   80
        100   60
        25   25
        75   25
        25   75
        60   82.5
        82.5  60];
    
    Elementos = [1  2  5  4 12;      %Matriz de conectividades Q5
        2  3  6  5 13;
        4  5  8  7 14;
        5 10  9  8 15;
        5  6 11 10 16];
    
elseif strcmp(eleType,'Q8') || strcmp(eleType,'AHMAD8')
    nNodEle = 8;
elseif strcmp(eleType,'Q9') || strcmp(eleType,'AHMAD9')
    nNodEle = 9;
end

NNodos = size(Nodos,1);
NEle = size(Elementos,1);

title('Malla')
meshplot(Elementos,Nodos,'b',1)

% Grados de Libertad
nDofNod = 2;
nDofTot = nDofNod*NNodos;

bc = zeros(NNodos,nDofNod);

bc(1,1:2) = 1;
bc([2,3],2) = 1;
bc = logical(reshape(bc',[],1));

% Cargas
R = zeros(NNodos,nDofNod);
R(7,2)= 2.5;
R(8,2)= 5;
R(9,2)= 2.5;

R = reshape(R',[],1);

%% Cálculos
rsInt = ones(1,2)*2; %*1,*2,*3
[wpg, upg, npg] = gauss(rsInt);

% Matriz de rigidez
K = zeros(nDofTot);
for iele = 1:NEle
    Ke = zeros(nNodEle*nDofNod);
    nodesEle = Nodos(Elementos(iele,:),:);   
    for ipg = 1:npg
        % Punto de Gauss
        ksi = upg(ipg,1);
        eta = upg(ipg,2);
        
        % Funciones de forma respecto de ksi, eta
        NF = shapefuns([ksi eta],eleType);
        N(1,1:2:nDofNod*nNodEle)=NF;
        N(2,2:2:nDofNod*nNodEle)=NF;
        
        % Derivadas de las funciones de forma respecto de ksi, eta
        dN = shapefunsder([ksi eta],eleType);
        % Derivadas de x,y, respecto de ksi, eta
        jac = dN*nodesEle(1:nNodEle,:);
        % Derivadas de las funciones de forma respecto de x,y.
        dNxy = jac\dN;          % dNxy = inv(jac)*dN
        
        B = zeros(size(C,2),nNodEle*nDofNod);
        B(1,1:2:nDofNod*nNodEle-1) = dNxy(1,:);
        B(2,2:2:nDofNod*nNodEle) = dNxy(2,:);
        B(3,1:2:nDofNod*nNodEle-1) = dNxy(2,:);
        B(3,2:2:nDofNod*nNodEle) = dNxy(1,:);
        
        Ke = Ke + t*B'*C*B*wpg(ipg)*det(jac);
    end
    eleDofs = node2dof(Elementos(iele,1:nNodEle),nDofNod);
    K(eleDofs,eleDofs) = K(eleDofs,eleDofs) + Ke;
end

% Reducción
Fixed = bc;
Free = ~bc;
Kr = K(Free,Free);
Rr = R(Free);

% Solve
Dr = Kr\Rr;

% Reconstrucción
D = zeros(length(bc),1);
D(Free) = Dr;

% Desplazamientos
Desp = reshape(D,nDofNod,NNodos)';
NodosDesp = Nodos+Desp;
NodosDespScale = Nodos+S*Desp;

%% Gráficos
title(['Deformada (escala: ',num2str(S),')'])
hold on
meshplot(Elementos,NodosDespScale,'k',0)

disp('Desplazamiento nodo 1: ')
Desp(1,:)

%% Recuperación de tensiones en los nodos
stress = zeros(NEle,nNodEle,3);
strain = zeros(NEle,nNodEle,3);
uNod = [ -1 -1
          1 -1
          1  1
         -1  1 ];

Ni = shapefuns   (uNod,eleType);
dN = shapefunsder(uNod,eleType);     
for iele = 1:NEle
    nodesEle = Nodos(Elementos(iele,:),:);
    for inode = 1:nNodEle
        % Derivadas de x,y, respecto de ksi, eta
        jac = dN(:,:,inode)*nodesEle;                      
        % Derivadas de las funciones de forma respecto de x,y.
        dNxy = jac\dN(:,:,inode);          % dNxy = inv(jac)*dN(:,:,ipg)
        r = Ni(:,:,inode)*nodesEle(:,1);
        
        B = zeros(size(C,2),nDofNod*nNodEle);
        B(1,1:2:7) = dNxy(1,:);
        B(2,2:2:8) = dNxy(2,:);
        B(3,1:2:7) = dNxy(2,:);
        B(3,2:2:8) = dNxy(1,:);
        
        eleDofs = node2dof(Elementos(iele,:),nDofNod);
        def = B*D(eleDofs);
        strain(iele,inode,:) = def;
        stress(iele,inode,:) = C*def;
    end
end
bandplot(Elementos,NodosDespScale,stress(:,:,3),[],'k');
