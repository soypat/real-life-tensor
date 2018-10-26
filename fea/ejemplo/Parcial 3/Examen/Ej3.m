%% Parcial 3
%% EJ 3
% Alumno: Andres Battisti

clear
clc
close all

%% Datos
t=1;%mm
E=210000; %MPa
nu=0.3;
% Modifique las propiedasdes del material
S = 100000; % escala

% eleType = 'Q5'; % Selecciono Q5

%Plane Strain
C = (E/((1+nu)*(1-2*nu)))*[1-nu    nu      0;
                             nu  1-nu      0;
                              0    0  0.5-nu];

%% Discretización
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

NNodos = size(Nodos,1);
NEle = size(Elementos,1);

title('Malla')
meshplot(Elementos,Nodos,'b',1)

% Grados de Libertad
nDofNod = 2;
nDofTot = nDofNod*NNodos;

bc = zeros(NNodos,nDofNod);

bc(1,1:2) = true;
bc([2,3],2) = true;
bc = logical(reshape(bc',[],1));

%% Cargas
R = zeros(NNodos,nDofNod);
% Modifico las cargas a las adecuadas
q = 5; %[N/mm]
R(10,1)= -q;

R = reshape(R',[],1);

%% Funciones de forma para Q5 (de shapefunsQ5.m)
dNfunc =@(x,y) [ (x*y^2)/2 + y/4 - 1/4, (x*y^2)/2 - y/4 + 1/4, (x*y^2)/2 + y/4 + 1/4, (x*y^2)/2 - y/4 - 1/4, -2*x*y^2
                 (y*x^2)/2 + x/4 - 1/4, (y*x^2)/2 - x/4 - 1/4, (y*x^2)/2 + x/4 + 1/4, (y*x^2)/2 - x/4 + 1/4, -2*x^2*y];
         
%% Cálculos
[wpg, upg, npg] = gauss([2 2]); % Integro con 2 puntos de gauss

% Matriz de rigidez
K = zeros(nDofTot);
storeKe = zeros(10,10,5);
for iele = 1:NEle
    Ke = zeros(nNodEle*nDofNod);
    nodesEle = Nodos(Elementos(iele,:),:);   
    for ipg = 1:npg
        % Punto de Gauss
        ksi = upg(ipg,1);
        eta = upg(ipg,2);
        % Derivadas de las funciones de forma respecto de ksi, eta
        dN = dNfunc(ksi,eta);
        % Derivadas de x,y, respecto de ksi, eta
        jac = dN*nodesEle;
        % Derivadas de las funciones de forma respecto de x,y.
        dNxy = jac\dN;          % dNxy = inv(jac)*dN
        
        B = zeros(size(C,2),nNodEle*nDofNod);
        B(1,1:2:nDofNod*nNodEle-1) = dNxy(1,:);
        B(2,2:2:nDofNod*nNodEle) = dNxy(2,:);
        B(3,1:2:nDofNod*nNodEle-1) = dNxy(2,:);
        B(3,2:2:nDofNod*nNodEle) = dNxy(1,:);
        
        Ke = Ke + t*B'*C*B*wpg(ipg)*det(jac);
    end
    storeKe(:,:,iele) = Ke;
    eleDofs = node2dof(Elementos(iele,1:nNodEle),nDofNod);
    K(eleDofs,eleDofs) = K(eleDofs,eleDofs) + Ke;
%     eig(Ke) % verifico los autovalores nulos 
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

disp(strcat(['Desplazamiento nodo B: ']))
fprintf('\n')
disp(strcat([num2str(Desp(5,:)),' mm en x e y respectivamente']))
