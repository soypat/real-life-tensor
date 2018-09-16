
% Patch Test - Q8 iosparamétrico.
clear
close all
format short g
clc
% UNIDADES N - mm
TestType = 'Sigma_x' ; % 'Sigma_y' ; %  ; %'Tau_xy'
Integracion = 'Subintegrada'; % 'Full'; %
StressCalc  = 'PG'; %'PG' ; %'Nodos'
Lados = 'Rectos'; % 'Curvos'; %

% Discretizacion
switch Lados
    case 'Rectos'
        nod = [ 0.0   0.0
            0.5   0.0
            1.0   0.0
            1.5   0.0
            2.0   0.0
            0.0   0.5
            1.125 0.375
            2.0   0.5
            0.0   1.0
            0.625 0.875
            1.25  0.75
            1.625 0.875
            2.0   1.0
            0.0   1.5
            1.125 1.375
            2.0   1.5
            0.0   2.0
            0.5   2.0
            1.0   2.0
            1.5   2.0
            2.0   2.0 ];          % Coordenadas nodales lados rectos
        
    case 'Curvos'
        nod = [ 0.0   0.0
            0.5   0.0
            1.0   0.0
            1.5   0.0
            2.0   0.0
            0.0   0.5
            1.0   0.375
            2.0   0.5
            0.0   1.0
            0.625 1.0
            1.25  0.75
            1.625 0.75
            2.0   1.0
            0.0   1.5
            1.0   1.375
            2.0   1.5
            0.0   2.0
            0.5   2.0
            1.0   2.0
            1.5   2.0
            2.0   2.0 ];          % Coordenadas nodales lados curvos
end

      elem = [ 1  3  11  9   2  7  10   6
               3  5  13  11  4  8  12  7
               9  11  19  17 10 15 18 14
               11  13  21  19  12  16  20  15];  %Matriz de conectividades: ojo el orden! muy importante!

nDofNod = 2;                    % Número de grados de libertad por nodo
nNodEl = 8;                     % Número de nodos por elemento
nElem = size(elem,1);           % Número de elementos
nNod = size(nod,1);             % Número de nodos
nDofTot = nDofNod*nNod;         % Número de grados de libertad


R = zeros(nNod,nDofNod);        % Vector de cargas
bc = false(nNod,nDofNod);       % Matriz de condiciones de borde
switch TestType 
    case 'Sigma_x' %sigma x constante
        
        bc(1,1:2) = true;
        bc([6 9 14],1) = true;
        
        
        R([5 21],1) = 1/6;
        R([8 16],1) = 2/3;
        R(13,1) = 1/3;
        R(17,1) = -1/6;
        
    case 'Sigma_y' %sigma y constante
      bc(1,1:2) = true;
        bc(2:4,2) = true;
        
        R(5,2) = -1/6;
        R([17 21],2) = 1/6;
        R([18 20],2) = 2/3;
        R(19,2) = 1/3; 
        
    case 'Tau_xy' %corte xy constante
        bc(1,1:2) = true;
        bc([6 9 14],1) = true;
        
        R([5 21],2)= 1/6;
        R([8 16],2)= 2/3;
        R(13,2) = 1/3;
        
        R(17,2)= -1/6;
        R([6 14],2) = -2/3;
        R(9,2) = -1/3;
        
        R([17 21],1) = 1/6;
        R([18 20],1) = 2/3;
        R(19,1) = 1/3;
        
        R(5,1)= -1/6;
        R([2 4],1)= -2/3;
        R(3,1) = -1/3;
        
end

meshplot(elem,nod,'b')
axis equal

% Propiedades del Material
E = 1;
NU = 0.33;
C = E/(1 - NU^2)*[ 1.0     NU         0.0
                    NU    1.0         0.0
                   0.0    0.0     (1 - NU)/2 ];

%% Matriz de rigidez

switch Integracion
    case 'Subintegrada'
        PG = [-1/sqrt(3) 1
             1/sqrt(3)   1]; %La segunda columna es el peso
    case 'Full'
        PG = [-sqrt(0.6) 5/9
             0           8/9
             sqrt(0.6)   5/9]; %La segunda columna es el peso
end

ordInt = size(PG,1); %Orden de integración Gaussiana

for ipg = 1:ordInt
    csi = PG(ipg,1);
    for jpg = 1:ordInt
        eta = PG(jpg,1);
        varName = genvarname(['dN' num2str(ipg) num2str(jpg)]); 
        %genvarname es una función de matlab que te permite guardar un
        %string como un nombre de variable. Así, puedo incluir los
        %contadores del for en el nombre de la variable que quiero guardar.
        %Guardo dN11, dN12, dN21, dN22
        eval([varName ' = dNQ8(csi, eta)']);
    end
end

%Guardé las derivadas de las funciones de forma respecto a csi y eta, que
%es la "parte independiente del Jacobiano", evaluada en los 4 o 9 puntos de
%integración Gaussiana. La llamo independiente porque va a ser igual para
%cualquier elemento, si bien depende de csi y eta, no depende en absoluto
%de x ni de y.

% Buscamos matriz de rigidez y ensamblamos
K = zeros(nDofTot);
nodeDofs = reshape(1:nDofTot,nDofNod,nNod)';
for iElem = 1:nElem
    %Matriz de rigidez elemental (elementos rectangulares)
    valNod = nod(elem(iElem,:),:);    
    ke = zeros(16);
    for ipg = 1:ordInt
        for jpg = 1:ordInt
        dN = eval(['dN' num2str(ipg) num2str(jpg)]); %El "Jacobiano independiente" evaluada en el PG que nos interesa
        jac = dN*valNod; %Jacobiano es el producto entre la "parte independiente" y los valores nodales.
        detJ = det(jac);
        dNxy = jac\dN; %Derivadas de las funciones de forma respecto a x e y.
        
        B = zeros(3,16);
        B(1,1:2:15) = dNxy(1,:);
        B(2,2:2:16) = dNxy(2,:);
        B(3,1:2:15) = dNxy(2,:);
        B(3,2:2:16) = dNxy(1,:); 
          
        w = [PG(ipg,2) PG(jpg,2)];
        ke = ke + w(1)*w(2)*B.'*C*B*detJ;
        end
    end
    eleDofs = reshape(nodeDofs(elem(iElem,:),:)',1,[]);
    K(eleDofs,eleDofs) = K(eleDofs,eleDofs) + ke;
end

%% Calculo desplazamientos

%Reducción matriz
isFixed = reshape(bc',[],1);
isFree = ~isFixed;

Rr = reshape(R',[],1);

%Solver
Dr = K(isFree,isFree)\Rr(isFree);

%Reconstrucción
D = zeros(nDofTot,1);
D(isFree) = Dr;

%Configuración deformada
Dreshape = reshape(D',2,[])';
meshplot(elem,nod + Dreshape,'r') %Dibujo deformación en rojo
axis equal

switch StressCalc
    case 'PG'
        %% Tensión en los puntos de superconvergencia
        
        PG = [-1/sqrt(3) 1
            1/sqrt(3)   1]; %buscamos un orden menor al ´full´, para lograr la superconvergencia
        
        %en el caso de integracion full, calculamos nuevamente el valor de los
        %jacobianos en los puntos de gauss, pq nos interesa la integracion
        %gausseana de un orden menor que el full
        %Extrapolamos multiplicando por raiz de 3 desde los puntos de
        %integracion gausseana hacia los nodos.
        if strcmp(Integracion, 'Full')
            for ipg = 1:2
                csi = PG(ipg,1);
                for jpg = 1:2
                    eta = PG(jpg,1);
                    varName = genvarname(['dN' num2str(ipg) num2str(jpg)]);
                    %genvarname es una función de matlab que te permite guardar un
                    %string como un nombre de variable. Así, puedo incluir los
                    %contadores del for en el nombre de la variable que quiero guardar.
                    %Guardo dN11, dN12, dN21, dN22
                    eval([varName ' = dNQ8(csi, eta)']);
                end
            end
        end
        
        % Buscamos matriz deformacion-desplazamientos
        Stress = zeros(8*nElem,3); %8 filas por cada elementos: las tres tensiones del nodo 1, las tres tensiones del nodo 2, etc.
        StressAvg = zeros(nNod,3);
        ind = zeros(nNod,1); %Cuenta cuantas veces está compartido un nodo, para después dividir por eso al promediar
        for iElem = 1:nElem
            valNod = nod(elem(iElem,:),:);
            PGstress = zeros(4,3); % 4 puntos PG, 3 tensiones
            eleDofs = reshape(nodeDofs(elem(iElem,:),:)',1,[]);
            % 4-------3
            % |       |
            % |       |
            % 1-------2 Estos son los PG de superconvergencia donde calculamos
            % tensiones. 11 es el PG 1, 12 es el PG 4, 21 es el PG 2, 22 es el PG
            % 3.
            RSextrapolation = [-sqrt(3)  -sqrt(3)
                sqrt(3)  -sqrt(3)
                sqrt(3)   sqrt(3)
                -sqrt(3)  sqrt(3)
                0      -sqrt(3)
                sqrt(3)     0
                0       sqrt(3)
                -sqrt(3)    0];
            for ipg = 1:2
                for jpg = 1:2
                    dN = eval(['dN' num2str(ipg) num2str(jpg)]); %El "Jacobiano independiente" evaluada en el PG que nos interesa
                    jac = dN*valNod; %Jacobiano es el producto entre la "parte independiente" y los valores nodales.
                    dNxy = jac\dN; %Derivadas de las funciones de forma respecto a x e y.
                    B = zeros(3,16);
                    B(1,1:2:15) = dNxy(1,:);
                    B(2,2:2:16) = dNxy(2,:);
                    B(3,1:2:15) = dNxy(2,:);
                    B(3,2:2:16) = dNxy(1,:);
                    
                    if ipg == jpg
                        PGstress(ipg + jpg - 1,:) = C*B*D(eleDofs);
                    else
                        PGstress(ipg * jpg^2,:) = C*B*D(eleDofs);
                    end
                    %Ahora que tenemos las tensiones en los puntos Gaussianos dentro
                    %del elemento, extrapolamos a los nodos y guardamos esa data
                    StressElem = zeros(8,3);
                    for iNod = 1:8
                        r = RSextrapolation(iNod,1);
                        s = RSextrapolation(iNod,2);
                        N1 = (1 - r)*(1 - s) / (4);
                        N2 = (1 + r)*(1 - s) / (4);
                        N3 = (1 + r)*(1 + s) / (4);
                        N4 = (1 - r)*(1 + s) / (4);
                        N = [N1 N2 N3 N4];
                        StressElem(iNod,:) = N * PGstress;
                    end
                end
            end
            Stress(8*iElem-7:8*iElem,:) = StressElem;
            % Las tensiones están ordenadas por elemento, queremos tenerlas en orden
            % por nodo, promediando entre los valores.
            NodosAfectados = elem(iElem,:);
            StressAvg(NodosAfectados,:) = StressAvg(NodosAfectados,:) + StressElem;
            ind(NodosAfectados,:) = ind(NodosAfectados,:) + 1;
        end
        for iNod = 1:nNod
            StressAvg(iNod,:) = StressAvg(iNod,:)./ind(iNod);
        end
    %% Tensión calculada en los nodos    
    case 'Nodos'
        % Para conocer las tensiones en los nodos, debemos evaluar la
        % matriz B en los nodos. Para eso, tenemos que conocer el valor del
        % jacobiano en los nodos.
        CsiEtaNod = [-1 -1
            1  -1
            1  1
            -1  1
            0  -1
            1  0
            0  1
            -1 0];
        Stress = zeros(8*nElem,3); %8 filas por cada elemento: las tres tensiones del nodo 1, las tres tensiones del nodo 2, etc.
        StressAvg = zeros(nNod,3);
        ind = zeros(nNod,1); %Cuenta cuantas veces está compartido un nodo, para después dividir por eso al promediar
        for i = 1:8
            csi = CsiEtaNod(i,1);
            eta = CsiEtaNod(i,2);
            varName = genvarname(['dN' num2str(i)]);
            %genvarname es una función de matlab que te permite guardar un
            %string como un nombre de variable. Así, puedo incluir los
            %contadores del for en el nombre de la variable que quiero guardar.
            %Guardo dN1, dN2, dN3, dN4, etc
            eval([varName ' = dNQ8(csi, eta)']);
        end
        for iElem = 1:nElem
            valNod = nod(elem(iElem,:),:);
            NodStress = zeros(8,3); % 8 nodos, 3 tensiones
            eleDofs = reshape(nodeDofs(elem(iElem,:),:)',1,[]);
            for i = 1:8
                dN = eval(['dN' num2str(i)]); %El "Jacobiano independiente" evaluada en el nodo que nos interesa
                jac = dN*valNod; %Jacobiano es el producto entre la "parte independiente" y los valores nodales.
                dNxy = jac\dN; %Derivadas de las funciones de forma respecto a x e y.
                B = zeros(3,16);
                B(1,1:2:15) = dNxy(1,:);
                B(2,2:2:16) = dNxy(2,:);
                B(3,1:2:15) = dNxy(2,:);
                B(3,2:2:16) = dNxy(1,:);
                
                NodStress(i,:) = C*B*D(eleDofs);
            end
            %  Ya tenemos las tensiones nodales del elemento. Ahora lo
            %  guardamos en la matriz global de tensiones.
            Stress(8*iElem-7:8*iElem,:) = NodStress;
            % Las tensiones están ordenadas por elemento, queremos tenerlas en orden
            % por nodo, promediando entre los valores.
            NodosAfectados = elem(iElem,:);
            StressAvg(NodosAfectados,:) = StressAvg(NodosAfectados,:) + NodStress;
            ind(NodosAfectados,:) = ind(NodosAfectados,:) + 1;
        end
        for iNod = 1:nNod
            StressAvg(iNod,:) = StressAvg(iNod,:)./ind(iNod);
        end
        
end

disp(StressAvg)

%% Ploteo de deformada
figure
for iElem = 1:nElem
    eleDofs = reshape(nodeDofs(elem(iElem,:),:)',1,[]); %Ubicacion de los dofs nodales
    dispElem = D(eleDofs);
    plotDef(nod(elem(iElem,:),:),dispElem,'r','Q8')
    hold on
end