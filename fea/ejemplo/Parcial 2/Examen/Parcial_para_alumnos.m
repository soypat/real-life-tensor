clear 
close all
clc

%% Datos Genéricos
t=10;    %mm
E=1000;    %MPa
nu=.45;
rho=1;  %kg/mm3

scale = 10;

%Plane Strain
% C = (E/((1+nu)*(1-2*nu)))*[1-nu    nu      0;
%                            nu  1-nu      0;
%                             0    0  0.5-nu];
%Plane Stress ya que tiene un espesor fino comparado con el di'ametro
C = (E/(1-nu^2))*[1 nu 0;
                  nu 1 0;
                  0 0 (1-nu)/2];

% Cargar malla
nNodEle = 8;
load('Malla2.mat');
Elementos = elem;
NEle = size(Elementos,1);
Nodos = nodes;
NNodos =size(Nodos,1);

figure(1)
title('Malla')
meshplot(Elementos,Nodos,'k')

% Grados de Libertad
nDofNod = 2;
nDofTot = nDofNod*NNodos;

horAxis=find(abs(Nodos(:,2))<1e-6);
verAxis=find(abs(Nodos(:,1))<1e-6);

bc = zeros(NNodos,nDofNod);
bc([9 21 6 20 3],2) = true;
bc([7 17 4 16 1],1) = true;
bc = logical(reshape(bc',[],1));

%% Cálculos
% Puntos de Gauss integración full
rsInt = ones(1,2)*3; %*1,*2,*3
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
        
        % Derivadas de las funciones de forma respecto de ksi, eta
        dN = shapefunsder([ksi eta],'Q8');  
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
                        
        % Funciones de forma evaluadas en coordenadas ksi y eta
        NF = shapefuns([ksi eta],'Q8');
    end
    eleDofs = node2dof(Elementos(iele,:),nDofNod);
    K(eleDofs,eleDofs) = K(eleDofs,eleDofs) + Ke;
end

%% Fuerzas de volumen
RG=zeros(NNodos,nDofNod);
omega = 100*pi()/30; %[rad/s]
rhopol = 900e-9; %[kg/m^3]
fx =@(x) rhopol*omega^2*x;
fy =@(y) rhopol*omega^2*y;
for iele = 1:NEle
    nodesEle = Nodos(Elementos(iele,:),:);
    for ipg = 1:npg
        ksi = upg(ipg,1);
        eta = upg(ipg,2);
        % Derivadas de las funciones de forma respecto de ksi, eta
        dN = shapefunsder([ksi eta],'Q8');  
        % Derivadas de x,y, respecto de ksi, eta
        jac = dN*nodesEle;
        N = shapefuns([ksi eta],'Q8');
        for n=1:8
            RG(Elementos(iele,:),1) = RG(Elementos(iele,:),1) + ...
                N(n)*fx(nodesEle(n,1))*t*wpg(ipg)*det(jac);
            RG(Elementos(iele,:),2) = RG(Elementos(iele,:),2) + ...
                N(n)*fy(nodesEle(n,2))*t*wpg(ipg)*det(jac);
        end
    end
end
RG = reshape(RG',[],1);
% Reducción
Fixed = bc;
Free = ~bc;
Kr = K(Free,Free);
Rr = RG(Free); %Elijo si volumen o presión	

    
% Solve
Dr = Kr\Rr;
    
% Reconstrucción Volumen
D = zeros(length(bc),1);
D(Free) = Dr;

% Desplazamientos
D = reshape(D,nDofNod,NNodos)';
NodosDesp = Nodos+D;
NodosDespScale = Nodos+scale*D;

%% Gráficos

figure(2)
title('Deformada')
hold on
meshplot(Elementos,Nodos,'k')
meshplotScale(Elementos,NodosDespScale,'r')
hold off