clear 
close all
clc
%% TEMA 1
%% Correcciones y explicación de mi código
% Mi idea fue crear un ciclo que vaya incrementando la el módulo de
% elasticidad del polímero, cree una matriz de rigidez y resuelva los
% desplazamientos hasta que se verifique que la posición final en
% cualquiera de los nodos que se encuentran en el diámetro interior llegue
% a 25mm, el punto donde no existirá más interferencia entre el acero y el
% polímero. Una vez encontrado este E, otro ciclo va incrementando la
% presión en el lugar de contacto con eje hasta que, nuevamente, el
% desplazamiento sea de 0.1mm.

% Me di cuenta tarde de que los nodos estaban mal recorridos, nunca me
% fijé si la matriz de rigidez estaba mal y, como los desplazamientos me
% daban algo razonable, supuse que las funciones de forma que dan el
% 'shapefuns' y shapefunsder' estaban en el mismo orden que la matriz
% Elementos. Ahora agregué una línea que corrige esa matriz.

% Corregí el vector de cargas. Yo había intentado cargarlo buscando la
% evaluación de las tensiones en los puntos de gauss cuando en realidad es
% más sencillo interpolar los mismos con las funciones de forma evaluadas
% en los puntos de gauss y las tensiones en los nodos (N*f(nodesEle))

% Terminé la parte del ciclo de tensiones superficiales que no llegué a
% completar.
%% Suposiciones
% Supongo que el desplazamiento del acero es mucho menor que el del
% polímero. Por lo tanto sólo calculo busco que la posición final de
% cualquiera de los nodos que se encuentran con el acero sea mayor que
% 25mm.

%ES LINEAL, NO HACIA FALTA CICLAR, BOLOOOOOOOO

%CON EL EJE INCLUIDO Y SUPONIENDO UN DELTA T SE PUEDE CALCULAR LA TENSION
% DE MONTAJE

%% Datos Genéricos

% Cargar malla
nNodEle = 8;
load('Malla2.mat');
Elementos = elem;
Elementos  = Elementos(:,[4 3 2 1 7 6 5 8]); % Correccción de elementos
NEle = size(Elementos,1);
Nodos = nodes;
NNodos =size(Nodos,1);

figure(1)
title('Malla')
meshplot(Elementos,Nodos,'k')

% Grados de Libertad
nDofNod = 2;
nDofTot = nDofNod*NNodos;

%horAxis=find(abs(Nodos(:,2))<1e-6);
%verAxis=find(abs(Nodos(:,1))<1e-6);

bc = zeros(NNodos,nDofNod);
bc([9 21 6 20 3],2) = true;
bc([7 17 4 16 1],1) = true;
bc = logical(reshape(bc',[],1));

%% Cálculos

for n = 1:199
% Puntos de Gauss integración full
rsInt = ones(1,2)*3; %*1,*2,*3
[wpg, upg, npg] = gauss(rsInt);

t=10;    %mm
E = 2000:-10:1;    %MPa
nu=.45;
scale = 100;

%Plane Stress ya que tiene un espesor fino comparado con el di'ametro
C = (E(n)/(1-nu^2))*[1 nu 0;
                  nu 1 0;
                  0 0 (1-nu)/2];

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
    centroelem = sum(nodesEle)/8; % Centro del elemento en coordenadas entructurales
        for ipg = 1:npg
            ksi = upg(ipg,1);
            eta = upg(ipg,2);
            x = centroelem(1)+ksi; % Coordenadas x e y de los puntos de gauss estructurales
            y = centroelem(2)+eta;
            % Derivadas de las funciones de forma respecto de ksi, eta
            dN = shapefunsder([ksi eta],'Q8');  
            % Derivadas de x,y, respecto de ksi, eta
            jac = dN*nodesEle;
            N = shapefuns([ksi eta],'Q8');
            Fx = N*fx(nodesEle(:,1));
            Fy = N*fx(nodesEle(:,2));
            RG(Elementos(iele,:),1) = RG(Elementos(iele,:),1) + ...
                N'*Fx*t*wpg(ipg)*det(jac);
            RG(Elementos(iele,:),2) = RG(Elementos(iele,:),2) + ...
                N'*Fy*t*wpg(ipg)*det(jac);
        end 
end
% Estaba evaluando mal las fuerzas

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

if NodosDesp(7,2) > 25;
    E = E(n)
    break
end
end

%% Gráficos

figure(2)
title('Deformada por fuerza centrífuga')
hold on
meshplot(Elementos,Nodos,'k')
meshplotScale(Elementos,NodosDespScale,'r')
hold off


% Para la tensi'on incial de montaje tengo que encontrar una presi'on en la
% superficie de contacto contra el eje que desplaze esos nodos hasta 25mm.
% Es decir, una presi'on de contacto equivalente a las fuerzas de volumen
% por centrifugaci'on. 

% Puntos de Gauss integración full
rsInt = ones(1,2)*3; %*1,*2,*3
[wpg, upg, npg] = gauss(rsInt);

t=10;    %mm
E = 2000:-10:1;    %MPa
nu=.45;
scale = 100;

%Plane Stress ya que tiene un espesor fino comparado con el di'ametro
C = (E(n)/(1-nu^2))*[1 nu 0;
                  nu 1 0;
                  0 0 (1-nu)/2];

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
    end
    eleDofs = node2dof(Elementos(iele,:),nDofNod);
    K(eleDofs,eleDofs) = K(eleDofs,eleDofs) + Ke;
end

%% Fuerzas de superficie
for n = 1:100
sig = 0:0.1:10; % [MPa] modificar hasta que de
sigma = sig(n);
RG=zeros(NNodos,nDofNod);
[wpg,upg] = gauss1D(2);
for iele = 1:NEle
    nodesEle = Nodos(Elementos(iele,:),:);
    if nodesEle(1,2) == 24.9 | nodesEle(2,1) == 24.9
        for ipg = 1:2
            ksi = upg(ipg);
            eta = -1;
            % Derivadas de las funciones de forma respecto de ksi, eta
            dN = shapefunsder([ksi eta],'Q8');  
            % Derivadas de x,y, respecto de ksi, eta
            jac = dN*nodesEle;
            N = shapefuns([ksi eta],'Q8');
            RG(Elementos(iele,1),1) = RG(Elementos(iele,1),1) + ...
                N(1)*sigma*t*wpg(ipg)*jac(1,1);
            RG(Elementos(iele,5),1) = RG(Elementos(iele,5),1) + ...
                N(5)*sigma*t*wpg(ipg)*jac(1,1);
            RG(Elementos(iele,2),1) = RG(Elementos(iele,2),1) + ...
                N(2)*sigma*t*wpg(ipg)*jac(1,1);
            RG(Elementos(iele,1),2) = RG(Elementos(iele,1),2) + ...
                N(1)*-sigma*t*wpg(ipg)*jac(1,2);
            RG(Elementos(iele,5),2) = RG(Elementos(iele,5),2) + ...
                N(5)*-sigma*t*wpg(ipg)*jac(1,2);
            RG(Elementos(iele,2),2) = RG(Elementos(iele,2),2) + ...
                N(2)*-sigma*t*wpg(ipg)*jac(1,2);
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

if NodosDesp(7,2) > 25;
    sigma
    break
end
end
%% Gráficos

figure(3)
title('Deformada por ensamble')
hold on
meshplot(Elementos,Nodos,'k')
meshplotScale(Elementos,NodosDespScale,'r')
hold off
