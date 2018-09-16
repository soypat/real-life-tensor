%% Correcciones
% Ortografía y aclaración del Ej 1
% Aclaraciones a lo largo de todo el código
% Borrado la matriz de DOF
% Modificado de orden de elementos
% Modificado de orden de grados de libertad, estaban muy desordenados
% Modificado: dofpornodo*nnod*3 --> Doftotal
% Modificado: se agregó la respuesta 2.c
%% EJ 1
% La matriz de rigidez debe ser simétrica ya que la relación entre los
% desplazamientos y las cargas es lineal, el material es linealmente elástico.
% Esto se puede demostrar aplicando  el teorema de Betti-Maxwell que se 
% encuentra en la página 37 del libro de Cook.
% La diagonal de la matriz debe ser positiva ya que los término Kii
% representan la relación directa entre los grados de libertad y las
% respectivas cargas que los afectan. Necesariamente el sentido en el
% que se desplazan o rotan los grados de libertad debe ser el mismo que las
% fuerzas o momentos aplicados sobre ellos.
%% EJ 2
clear; close all; clc
%% Datos
a = 40;          %mm
E = 2E5;         %MPa
L = 1000;        %mm
F = 1000;        %N
d = 40;          %mm
Av = a^2;        %mm^2
Ap = d^2*pi()/4; %mm^2
A = [Av Ap Av Av Av Av]; % Vector de areas por elemento
Iv = a^4/12;             % Momento de inercia en z para vigas
Ip = d^4*pi()/64;        %Momento de inercia en z para el pistón, no necesario

%% Nodos y elementos
cnod = [0 0
        sqrt(2)/4 sqrt(2)/4
        0 sqrt(2)/2
        sqrt(2)/2 sqrt(2)/2
        1+sqrt(2)/2 sqrt(2)/2
        2+sqrt(2)/2 sqrt(2)/2
        2 0]*L; % Coordenadas de los nodos
nnod = 7;       % Número de nodos
elem = [1 2
        2 3
        2 4
        4 5
        5 6
        6 7]; % Matriz de elementos
% Los elementos [1 3 4 5] son vigas ya que transmite momentos en sus
% extremos, los elementos [2 6] son barras.
nelem = 6;      % Número de elementos

%% DOF
dofpornodo = 3; % Grados de libertad por nodo, los lugares reticulados queda libre la rotación.
reticulas = 3; % Cantidad de puntos entre elementos que estan reticulados
Doftotal = dofpornodo*nnod+reticulas; % Grados de libertad totales

%% Matriz de rigidez
% Los elementos barras fueron conformados con la misma matriz de rigidez
% que las vigas. Esto crea grados de libertad inecesarios que perjudican la
% velocidad del programa pero en el momento del examen no di a tiempo para
% preparar uno separado.
% La matriz fue armada manualmente, elemento a elemento.
% Se agregó un grado de libertad por cada unión reticulada entre elementos,
% un total de 3. En el programa DOFs se puede ver la asignación
% numérica que tiene cada grado de libertad, el esquema de la figura DOFs
% muestra donde se encuntra cada uno.

Kglobal = zeros(Doftotal);
% Elemento 1, viga
Kglobal([1 2 3 4 5 6],[1 2 3 4 5 6]) = Kglobal([1 2 3 4 5 6],[1 2 3 4 5 6])+viga(A(1),Iv,E,cnod,elem,1);
% Elemento 2, barra (pistón)
Kglobal([4 5 7 8 9 10],[4 5 7 8 9 10]) = Kglobal([4 5 7 8 9 10],[4 5 7 8 9 10])+viga(A(2),Ip,E,cnod,elem,2);
% Elemento 3, viga
Kglobal([4 5 6 11 12 13],[4 5 6 11 12 13]) = Kglobal([4 5 6 11 12 13],[4 5 6 11 12 13])+viga(A(3),Iv,E,cnod,elem,3);
% Elemento 4, viga
Kglobal([11 12 14 15 16 17],[11 12 14 15 16 17]) = Kglobal([11 12 14 15 16 17],[11 12 14 15 16 17])+viga(A(4),Iv,E,cnod,elem,4);
% Elemento 5, viga
Kglobal([15 16 17 18 19 20],[15 16 17 18 19 20]) = Kglobal([15 16 17 18 19 20],[15 16 17 18 19 20])+viga(A(5),Iv,E,cnod,elem,5);
% Elemento 6, barra
Kglobal([18 19 21 22 23 24],[18 19 21 22 23 24]) = Kglobal([18 19 21 22 23 24],[18 19 21 22 23 24])+viga(A(6),Iv,E,cnod,elem,6);

%% RESPUESTA a.
Kglobal;

%% BC
fijo = zeros(Doftotal,1);           % Vector de restringidos
fijo([1 2 8 9 22 23]) = ones(6,1);  % DOFs restringidos a desp o rot 0
libre = ~fijo;                      % Dofs con movimiento permitido

%% Cargas
P = zeros(Doftotal,1);              % Vector de carga
P(16) = -F;                         % Carga puntual

%% Reducción
Kglobalred = Kglobal(libre,libre);
Pred = P(libre);

%% Solver
Dred = Kglobalred\Pred;
D = zeros(Doftotal,1);
D(libre) = Dred;                    % Vector de desplazamientos completo

%% Tensiones en la sección B
Delem4 = D([11 12 14 15 16 17]);    % Desplazamientos del elemento 4
Bvigas = [ - 6/(L^2) + 12*L/2/(L^3)
                - 4/L + 6*L/2/(L^2)
             6/(L^2) - 12*L/2/(L^3)
                 -2/L + 6*L/2/(L^2)]';

SigmaFlexSup = a/2*E*Bvigas*Delem4([2 3 5 6]);
SigmaFlexInf = -SigmaFlexSup;
Bbarras = [-1/L 1/L];
SigmaAxial = E*Bbarras*Delem4([1 4]);

%% RESPUESTA b.
SigmaTotal = [SigmaFlexSup+SigmaAxial SigmaAxial SigmaFlexInf+SigmaAxial]';
% Tensiones en la parte superior, centro e inferior.

%% RESPUESTA c.
Reaccion4 = viga(A(4),Iv,E,cnod,elem,4)*D([11 12 14 15 16 17]);
FuerzaA = norm(Reaccion4([1 2]));
% Fuerza de corte en el perno A
%% RESPUESTAS
Kglobal;
SigmaTotal; %MPa
FuerzaA;    %N