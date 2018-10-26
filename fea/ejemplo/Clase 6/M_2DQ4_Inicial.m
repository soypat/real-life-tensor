% Modelo Q4 Tracción simple.
clear
close all
format short g
clc

syms x y as bs E NU real;

% Discretizacion
nodes = [ 0.0  0.0
          1.5  0.0
          3.0  0.0
          0.0  1.0
          1.5  1.0
          3.0  1.0
          0.0  2.0  
          1.5  2.0
          3.0  2.0 ];          % Coordenadas nodales

elements = [ 1  2  5  4
             2  3  6  5 
             4  5  8  7
             5  6  9  8 ];      %Matriz de conectividades

nDofNod = 2;                    % Número de grados de libertad por nodo
nNodEl = 4;                     % Número de nodos por elemento
nel = size(elements,1);         % Número de elementos
nNod = size(nodes,1);           % Número de nodos
nDofTot = nDofNod*nNod;         % Número de grados de libertad
t=0.5;                          % Espesor
bc = false(nNod,nDofNod);       % Matriz de condiciones de borde
bc(1,1:2) = true;
bc(4,1:2) = true;
bc(7,1:2) = true;

bc(2,2) = true;
bc(3,2) = true;

R = zeros(nNod,nDofNod);        % Vector de cargas
R(7,2) = -225*t;
R(8,2) = -225*t - 225*t;
R(9,2) = -225*t - 1000;

figure(5)
meshplot(elements,nodes,'b')

% Propiedades del Material
E = 30e6;
NU = 0.25;
% Plane stress
C = t*E/(1 - NU^2)*[ 1.0     NU         0.0
                      NU     1.0        0.0
                     0.0     0.0     (1 - NU)/2 ];

% Funciones de forma
N1 = (as - x)*(bs - y) / (4*as*bs);
N2 = (as + x)*(bs - y) / (4*as*bs);
N3 = (as + x)*(bs + y) / (4*as*bs);
N4 = (as - x)*(bs + y) / (4*as*bs);
N = [N1 N2 N3 N4];
% Derivadas de las funciones de forma respecto de x
dNx = diff(N,x);
% Derivadas de las funciones de forma respecto de y
dNy = diff(N,y);
% Matriz B simbólica
Bs = [ dNx(1)     0   dNx(2)      0  dNx(3)      0  dNx(4)      0
           0  dNy(1)      0   dNy(2)     0   dNy(3)     0   dNy(4)
       dNy(1) dNx(1)  dNy(2)  dNx(2) dNy(3)  dNx(3) dNy(4)  dNx(4) ];

K = zeros(nDofTot);
nodeDofs = reshape(1:nDofTot,nDofNod,nNod)';
for iele = 1:nel
    % Matriz de rigidez elemental (elementos rectangulares)
    a = abs((nodes(elements(iele,2),1) - nodes(elements(iele,1),1))) / 2;
    b = abs((nodes(elements(iele,3),2) - nodes(elements(iele,2),2))) / 2;
    B = subs(Bs,{as,bs},{a,b});
    Ke = subs(int(int(B'*C*B,x,-a,a),y,-b,b));
    eleDofs = nodeDofs(elements(iele,:),:);
    eleDofs = reshape(eleDofs',[],1);
    K(eleDofs,eleDofs) = K(eleDofs,eleDofs) + Ke;  
end

% Reduccion Matriz
isFixed = reshape(bc',[],1);
isFree = ~isFixed;

Rr = reshape(R',[],1);

% Solver
Dr = K(isFree,isFree)\Rr(isFree);

% Reconstrucción
D = zeros(nDofTot,1);
D(isFree) = D(isFree) + Dr;

% Reacciones
Rv = K(isFixed,isFree)*D(isFree);
reacciones=zeros(nDofTot,1);
reacciones(isFixed) = Rv;
reacciones = (reshape(reacciones,nDofNod,[]))';

%Recuperación de tensiones
stress = zeros(nel,nNodEl,3);
for iele = 1:nel
    figure(iele)
    a = abs((nodes(elements(iele,2),1) - nodes(elements(iele,1),1)) / 2);
    b = abs((nodes(elements(iele,3),2) - nodes(elements(iele,2),2)) / 2);
    B = subs(Bs,{as,bs},{a,b});
    eleDofs = nodeDofs(elements(iele,:),:);
    eleDofs = reshape(eleDofs',[],1);
    e = B*round(D(eleDofs)*1e3)/1e3;
    stressFun = C*e;
    ezcontourf(char(eval(stressFun(1))),[-a a -b b])
    axis equal
end

% Configuración deformada
D = (reshape(D,nDofNod,[]))';
nodePosition = nodes + D(:,1:2);

% %Graficación
% bandplot(elements,nodePosition,stress(:,:,2),[-1 1],'k');
figure(5)
hold on
meshplot(elements,nodePosition,'k')

