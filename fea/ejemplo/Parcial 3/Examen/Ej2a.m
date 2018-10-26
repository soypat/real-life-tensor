%% Parcial 3
%% EJ 2a
% Alumno: Andres Battisti

clear
clc
close all

%% Datos
E=210000; %MPa
nu=0.3;

S = 100000; % escala

eleType = 'Q4'; % Selecciono Q4

%% Constitutiva para axisimetricos

f = nu/(1-nu);
g = (1-2*nu)/(2*(1-nu));
C = [1 f f 0
     f 1 f 0
     f f 1 0
     0 0 0 g]*(1-nu)*E/((1+nu)*(1-2*nu));

%% Discretización
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

NNodos = size(Nodos,1);
NEle = size(Elementos,1);

title('Malla')
meshplot(Elementos,Nodos,'b',1)

% Grados de Libertad
nDofNod = 2;
nDofTot = nDofNod*NNodos;

%% Condiciones de borde
bc = zeros(NNodos,nDofNod);
% El eje de simetr'ia se encuentra restingido en r (x) y la base en z (y)
bc(1:3,2) = true;
bc([1 4 7],1) = true;
bc = logical(reshape(bc',[],1));


%% Cargas
R = zeros(NNodos,nDofNod);
% Modifico las cargas a las adecuadas
q = 5; %[N/mm]
rcarga = 100; % [mm] del eje de simentr'ia
R(6,1)= -q*2*rcarga*pi;
% La carga se encuentra distribuida alrededor de la barra

R = reshape(R',[],1);

%% Cálculos
[wpg, upg, npg] = gauss([3 3]);

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
        N = shapefuns([ksi eta],eleType);
        % Derivadas de las funciones de forma respecto de ksi, eta
        dN = shapefunsder([ksi eta],eleType);
        % Derivadas de x,y, respecto de ksi, eta
        jac = dN*nodesEle(1:nNodEle,:);
        % Derivadas de las funciones de forma respecto de x,y.
        dNxy = jac\dN;          % dNxy = inv(jac)*dN
        % Radio en los puntos evaluados
        r = N*nodesEle(:,1);
        
        B = zeros(size(C,2),nNodEle*nDofNod);
        B(1,1:2:nDofNod*nNodEle-1) = dNxy(1,:);
        B(2,1:2:nDofNod*nNodEle-1) = N/r;
        B(3,2:2:nDofNod*nNodEle)   = dNxy(2,:);
        B(4,1:2:nDofNod*nNodEle-1) = dNxy(2,:);
        B(4,2:2:nDofNod*nNodEle)   = dNxy(1,:);
        
        Ke = Ke + B'*C*B*wpg(ipg)*det(jac)*r*2*pi;
    end
    eleDofs = node2dof(Elementos(iele,1:nNodEle),nDofNod);
    K(eleDofs,eleDofs) = K(eleDofs,eleDofs) + Ke;
%     eig(Ke)
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


%% Recuperación de tensiones en los nodos
stress = zeros(NEle,nNodEle,4);
strain = zeros(NEle,nNodEle,4);
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
        % Radio en los puntos evaluados
        N = Ni(:,:,inode);
        r = N*nodesEle(:,1);
        
        B = zeros(size(C,2),nNodEle*nDofNod);
        B(1,1:2:nDofNod*nNodEle-1) = dNxy(1,:);
        B(2,1:2:nDofNod*nNodEle-1) = N/r;
        B(3,2:2:nDofNod*nNodEle)   = dNxy(2,:);
        B(4,1:2:nDofNod*nNodEle-1) = dNxy(2,:);
        B(4,2:2:nDofNod*nNodEle)   = dNxy(1,:);
        
        eleDofs = node2dof(Elementos(iele,:),nDofNod);
        def = B*D(eleDofs);
        strain(iele,inode,:) = def;
        stress(iele,inode,:) = C*def;
    end
end

%% Recuperación de tensiones en el centro de cada elemento
stressCentro = zeros(NEle,nNodEle,4);
strainCentro = zeros(NEle,nNodEle,4);
uNod = [0 0];

Ni = shapefuns   (uNod,eleType);
dN = shapefunsder(uNod,eleType);     
for iele = 1:NEle
    nodesEle = Nodos(Elementos(iele,:),:);
    for inode = 1:nNodEle
        % Derivadas de x,y, respecto de ksi, eta
        jac = dN*nodesEle;                      
        % Derivadas de las funciones de forma respecto de x,y.
        dNxy = jac\dN;          % dNxy = inv(jac)*dN(:,:,ipg)
        % Radio en los puntos evaluados
        N = Ni;
        r = N*nodesEle(:,1);
        
        B = zeros(size(C,2),nNodEle*nDofNod);
        B(1,1:2:nDofNod*nNodEle-1) = dNxy(1,:);
        B(2,1:2:nDofNod*nNodEle-1) = N/r;
        B(3,2:2:nDofNod*nNodEle)   = dNxy(2,:);
        B(4,1:2:nDofNod*nNodEle-1) = dNxy(2,:);
        B(4,2:2:nDofNod*nNodEle)   = dNxy(1,:);
        
        eleDofs = node2dof(Elementos(iele,:),nDofNod);
        def = B*D(eleDofs);
        strainCentro(iele,inode,:) = def;
        stressCentro(iele,inode,:) = C*def;
    end
end

%% Promediado de tensiones en los nodos
stressAvg = zeros(NNodos,4);
for inode = 1:NNodos
    [i,j] = find(Elementos == inode);
    nShare = length(i);
    for ishare = 1:nShare
        stressAvg(inode,:) = stressAvg(inode,:) + ...
            squeeze(stressCentro(i(ishare),j(ishare),:))';
    end
    stressAvg(inode,:) = stressAvg(inode,:)/nShare;
end

%% Tensiones en el elemento 1
% Como el primer elemento esta en contacto con el eje de simetr'ia no se
% pueden evaluar las tensiones en los nodos que conforman dicho eje. Por lo
% tanto, se pueden calcular las tensiones en los nodos evaluando en el
% centro de cada elemento (subintegrado) cada integral. Luego se pueden
% promediar las tensiones de cada nodo seg'un la nformaci'on de cada
% elemento. Lamentablemente, los nodos que se encuentran en el eje no
% tienen ning'un vecino para ser promediados.
disp('Tensiones promediadas S_r S_tita S_z T_zr:')
disp([ [1 2 5 4]' stressAvg([1 2 5 4],:)])
%% Bandplots
figure
bandplot(Elementos,NodosDespScale,stressCentro(:,:,2),[],'k');
figure
Elemento1 = [1 2 3 4]';
Nodos1 = Nodos(Elementos(1,:),:);
% bandplot(Elemento1,Nodos1,stressAvg(Elementos(1,:),2),[],'k');