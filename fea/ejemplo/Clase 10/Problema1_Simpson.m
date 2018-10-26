clear 
close all
clc


%% Datos

H = 1;       % Altura en mm
L = 1;       % Base en mm
t = 1;       % Ancho en mm
E = 5e9;     % en MPa
NU = 0.3;

C = (E/((1-NU*2)*(1+NU)))*[ 1-NU NU 0 ; NU 1-NU 0 ; 0 0 (1 - 2*NU)/2 ]; %Plane Strain

lambda = E*NU / ((1+NU)*(1-2*NU));
MU = E / (2*(1+NU));

f1 = 2*(MU+lambda) * 1e-2;
f2 = -(MU+lambda) * 1e-2;
g11 = @(y) ((2*MU + lambda)*y - 2*lambda)*1e-2;
g21 = @(y) MU*(1-2*y)*1e-2;
g12 = @(x) MU*(x-2)*1e-2;
g22 = @(x) (-2*(2*MU+lambda)*x + lambda) *1e-2;

FAD = 2; %Factor de Amplificación de Desplazamientos: Evidencia desplazamientos en los gráficos

% Malla
NroEleL = 1;
NroEleH = 1;
NroEle = NroEleL*NroEleH;

% Elementos
NroNodosEle = 4;
Elementos = EleMalladorQ4(NroEleL,NroEleH,NroEle);
AEle = (H*L)/NroEle;

% Nodos
Nodos = NodMalladorQ4(L,H,NroEleL,NroEleH);
NroNodos = max(max(Elementos));

% Grados de Libertad
GdL = 2;
GdLTot = GdL*NroNodos;
GdLNodos=reshape(1:GdLTot,GdL,NroNodos)';
bc = zeros(NroNodos,GdL);
for iNod=1:NroNodos
    if abs(Nodos(iNod,1))<1e-9
        bc(iNod,1:2) = 1;
    end
    if abs(Nodos(iNod,2))<1e-9
        bc(iNod,1:2) = 1;
    end
end
bc = logical(reshape(bc',[],1)); 

% Cargas(SIMPSON)
% Puntos de Gauss
a = 1/sqrt(3);
UbiPG = [ -a -a ; a -a ; a a ; -a a ];
NroPG = size(UbiPG,1);
WeiPg = ones(NroPG,1);

% Fuerzas de Superficie (SIMPSON)
R = zeros(NroNodos,GdL);
for iEle = 1:NroEle
    NodEle = Elementos(iEle,:);
    NodEle(5) = NodEle(1); 
    for iNod = 1:4
        
        % Lado Derecho
        if abs(Nodos(NodEle(iNod),1)-L)<1e-9
            if abs(Nodos(NodEle(iNod+1),1)-L)<1e-9   % Si cumple ambas es lado derecho
                nodA = NodEle(iNod);
                nodB = NodEle(iNod + 1);
                
                if Nodos(nodA,2)<Nodos(nodB,2)   % Me fijo quien esta primero para integrar
                    a = Nodos(nodA,2);
                    b = Nodos(nodB,2);
                else
                    a = Nodos(nodB,2);
                    b = Nodos(nodA,2);
                    k = nodB;
                    nodB = nodA;
                    nodA = k;
                end
                
                l = b-a;
                R(nodA,1) = l/6 * (g11(a) + 2*g11((a+b)/2)) + R(nodA,1);    % Regla de simpson
                R(nodA,2) = l/6 * (g21(a) + 2*g21((a+b)/2)) + R(nodA,2); 
                R(nodB,2) = l/6 * (g21(b) + 2*g21((a+b)/2)) + R(nodB,2); 
                R(nodB,1) = l/6 * (g11(b) + 2*g11((a+b)/2)) + R(nodB,1);              
            end
        end
        
        % Lado Superior
        if abs(Nodos(NodEle(iNod),2)-H)<1e-9
            if abs(Nodos(NodEle(iNod+1),2)-H)<1e-9    % Si cumple ambas es lado superior
                nodA = NodEle(iNod);
                nodB = NodEle(iNod + 1);
                
                if Nodos(nodA,1)<Nodos(nodB,1)   % Me fijo quien esta primero para integrar
                    a = Nodos(nodA,1);
                    b = Nodos(nodB,1);
                else
                    a = Nodos(nodB,1);                    
                    b = Nodos(nodA,1);
                    k = nodB;
                    nodB = nodA;
                    nodA = k;
                end
                
                l = b-a;
                R(nodA,1) = l/6 * (g12(a) + 2*g12((a+b)/2)) + R(nodA,1);    % Regla de simpson
                R(nodA,2) = l/6 * (g22(a) + 2*g22((a+b)/2)) + R(nodA,2); 
                R(nodB,2) = l/6 * (g22(b) + 2*g22((a+b)/2)) + R(nodB,2); 
                R(nodB,1) = l/6 * (g12(b) + 2*g12((a+b)/2)) + R(nodB,1);              
            end
        end
    end
end

% Fuerzas Volumétricas
for iEle=1:NroEle
    for iNod=1:NroNodosEle
        R(Elementos(iEle,iNod),:) = R(Elementos(iEle,iNod),:)+[f1*AEle/4 f2*AEle/4];
    end
end    
R = reshape(R',[],1);


%% Cálculos

% Matriz de Rigidez
K = zeros(GdLTot);
for iEle=1:NroEle
    Ke = zeros(GdL*NroNodosEle);
    NodEle = Nodos(Elementos(iEle,:),:);
    for iPG=1:NroPG
        ksi = UbiPG(iPG,1);
        eta = UbiPG(iPG,2);
        dN = 1/4*[-(1-eta)   1-eta   1+eta -(1+eta)
                  -(1-ksi) -(1+ksi)  1+ksi   1-ksi ];
        Jac = dN*NodEle;
        dNxy = Jac\dN;
        
        B = zeros(size(C,2),GdL*NroNodosEle);
        B(1,1:2:7) = dNxy(1,:);
        B(2,2:2:8) = dNxy(2,:);
        B(3,1:2:7) = dNxy(2,:);
        B(3,2:2:8) = dNxy(1,:);
        
        Ke = Ke + B'*C*B*WeiPg(iPG)*det(Jac)*t;
    end
    GdLEle = reshape((GdLNodos(Elementos(iEle,:),:))',[],1);
    K(GdLEle,GdLEle) = K(GdLEle,GdLEle) + Ke;
end

% Reducción
Fijo = bc;
Libre = ~bc;
Kr = K(Libre,Libre);
Rr = R(Libre);	
    
% Solve
Dr = Kr\Rr;
    
% Reconstrucción
De = zeros(length(bc),1);
De(Libre) = De(Libre)+Dr;
    
% Desplazamientos
D = reshape(De,GdL,NroNodos)';
CD = Nodos+D(1:NroNodos,1:2);
CDG = Nodos+FAD*D(1:NroNodos,1:2);

% Chequeo Desp
Ux= @(x,y) x*y * 1e-2;
Uy= @(x,y) -2*x*y * 1e-2;

Ucorr = zeros(NroNodos,GdL); 
for iNod = 1:NroNodos
    xnod = Nodos(iNod,1);
    ynod = Nodos(iNod,2);
    Ucorr(iNod,1) = Ux(xnod,ynod);
    Ucorr(iNod,2) = Uy(xnod,ynod);
end


%% Gráficos

figure('Name','Resultado Simpson','NumberTitle','off')
hold on
meshplot(Elementos,Nodos,'b')
meshplot(Elementos,CDG,'r')
hold off
title('Resultado Simpson')

figure('Name','Resultado Exacto','NumberTitle','off')
hold on
meshplot(Elementos,Nodos,'b')
meshplot(Elementos,Nodos+Ucorr,'r')
hold off
title('Resultado Exacto')

