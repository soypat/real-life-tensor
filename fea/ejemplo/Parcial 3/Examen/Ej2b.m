%% Parcial 3
%% EJ 2b
% Alumno: Andres Battisti

clear
clc
close all

%% Datos
E=210000; %MPa
nu=0.3;

S = 1000000; % escala

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
p = 1; %[N/mm]
rcarga = 80; % [mm] del eje de simentr'ia
R(11,:)= -p*2*rcarga*pi*sqrt(2)/2;
R(10,:)= -p*2*rcarga*pi*sqrt(2);
R(9,:)= -p*2*rcarga*pi*sqrt(2)/2;
% La carga se encuentra distribuida alrededor de la barra

R = reshape(R',[],1);

%% Cálculos
[wpg, upg, npg] = gauss([2 2]);

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

%% Rigid links
nodeDofs = reshape(1:nDofTot,nDofNod,NNodos)';

c = zeros(3,nDofTot);
u11 = nodeDofs(11,1);
u10 = nodeDofs(10,1);
u9 = nodeDofs(9,1);
v11 = nodeDofs(11,2);
v10 = nodeDofs(10,2);
v9 = nodeDofs(9,2);

c(1,[u9 u11 u10]) = [1 1 -2];
c(2,[v10 v11 u11 u10]) = [1 -1 1 -1];
c(3,[v9 v11 u11 u10]) = [1 -1 2 -2];
% c(4,[u10 u11]) = [1 -1];
% c(5,[v11 u11]) = [1 -1];
% c(1,[v11 u11]) = [1 -1];
% c(2,[v10 u11]) = [1 -1];
% c(3,[v9 u11]) = [1 -1];
% c(4,[u10 u11]) = [1 -1];
% c(5,[u9 u11]) = [1 -1];

Q = zeros(3,1);

constrained = false(nDofTot,1);
constrained([u9 v11 v9]) = true;
nconst = length(find(constrained));
Kc = [K c'
      c zeros(nconst)];
Rc = [ R; Q];


%% Reducción
Fixed = bc;
Free = ~bc;
Free = [ Free; true(nconst,1)];
Kr = Kc(Free,Free);
Rr = Rc(Free);

% Solve
Dr = Kr\Rr;

% Reconstrucción
D = zeros(length(bc),1);
D(Free) = Dr;

% Desplazamientos
Desp = reshape(D,nDofNod,NNodos)';
NodosDesp = Nodos+Desp;
NodosDespScale = Nodos+S*Desp;

%% Soluci'on
% No encuentro donde esta el error, seguramente es algebraico en las
% condiciones que impuse. Lo que estoy intentando hacer es que los 3 nodos
% 11 10 y 9 se esten vinculados rigidamente. Una vez calculado el
% desplazamiento mi intenci'on era aprovechar la linealidad del problema
% para calcular cuanto debe ser la carga distribuida 'p' para que se
% dezplace 2mm en la direcci'on impuesta a 45`.

%% Gráficos
title(['Deformada (escala: ',num2str(S),')'])
hold on
meshplot(Elementos,NodosDespScale,'k',0)
