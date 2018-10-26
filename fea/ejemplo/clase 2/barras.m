% Modelo de barras plano.
clear
clc
close all

% % Tranquera
% % Discretizacion
% nodes = [ 0.0     0.0     
%           1.5     0.0 
%           1.5     2.5
%           0.0     2.5 ];
%      
% elements = [ 1     2     
%              2     3 
%              3     4
%              4     1
%              1     3
%              2     4 ];         %Conectividades de elementos
%          
% A = [100 100 100 100 100 100];  % Area de cada elemento 
% 
% 
% nDofNod = 2;                    % Número de grados de libertad por nodo
% nel = size(elements,1);         % Número de elementos
% nNod = size(nodes,1);           % Número de nodos
% 
% bc = false(nNod,nDofNod);       % Matriz de condiciones de borde
% bc(1,1:2) = true;
% bc(2,2) = true;
% 
% R = zeros(nNod,nDofNod);               % Vector de cargas                   
% R([3 4],2) = 100000000;    % tracción simple

% Escuadra
% Discretizacion
L = 300;
alpha = 30;
P = 400000;
nodes = [ 0.0     0.0     
            L     0.0 
          0.0     L*tand(alpha) ];
     
elements = [ 2 3
             1 2 ];         %Conectividades de elementos
         
A = [5 10];  % Area de cada elemento 


nDofNod = 2;                    % Número de grados de libertad por nodo
nel = size(elements,1);         % Número de elementos
nNod = size(nodes,1);           % Número de nodos

bc = false(nNod,nDofNod);       % Matriz de condiciones de borde
bc(1,1:2) = true;
bc(3,1:2) = true;

R = zeros(nNod,nDofNod);               % Vector de cargas                   
R(2,2) = -P;

figure; hold on;
Draw_Barra(elements,nodes,'b');

disp('Configuración inicial')
disp(nodes)
disp('Conectividades')
disp(elements)
disp('Condiciones de Borde')
disp(bc)
disp('Cargas')
disp(R)

% Propiedades del Material
E = 2E6;

% Armado Matriz de Rigidez
K = zeros(nDofNod*nNod);

for iele = 1:nel;
    dir = nodes(elements(iele,2),:) - nodes(elements(iele,1),:);
    le = norm(dir);
    dir = dir/le;
    T = [ dir 0 0
          0 0 dir ];
    Ke =  A(iele)*E/le * [ 1 -1
                          -1  1 ];
    Ke = T'*Ke*T;
    eleDofs = node2dof(elements(iele,:),nDofNod);     
    K(eleDofs,eleDofs) = K(eleDofs,eleDofs) + Ke;
end

% Reduccion Matriz
fixed = reshape(bc',[],1);
free = ~fixed;

Rr = reshape(R',[],1);

% Solver
Dr = K(free,free)\Rr(free);

% Reconstruccion
D = zeros(nDofNod*nNod,1);
D(free) = D(free) + Dr;

% Verificar reacciones

% Tensiones
S = zeros(1,nel);
for iele = 1:nel;
    dir = nodes(elements(iele,2),:) - nodes(elements(iele,1),:);
    le = norm(dir);
    dir = dir/le;
    T = [ dir 0 0
          0 0 dir ];
    B = [-1 1]/le;
    eleDofs = node2dof(elements(iele,:),nDofNod);
    S(iele) = E * B * T*D(eleDofs);
end

defNodes = nodes + (reshape(D,nDofNod,[]))';

% Salida de datos
disp('Desplazamientos')
disp(D)

disp('Configuración deformada')
disp(defNodes)

disp('Tensiones normales')
disp(S')

Draw_Barra(elements,defNodes,'k')




