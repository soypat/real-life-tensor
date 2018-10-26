%% EJ 1
% Debe ser singular ya que no posee la información de las condiciónes de
% borde en la cual se encuentra el modelo a la cual pertenece. Las
% condiciones de borde son las que limitan al sistema de ecuaciones a tener
% una única solución. De no existir condiciones de borde o ser erroneas,
% por ejemplo un mecanismo, existen infinitas soluciones para el sistema de
% ecuaciones. Una matriz que respresenta un sistema de infinitas soluciones
% no es inversible.
%% EJ 2
% Un modelo físico es la forma de representar la realidad del problema, con
% todo detalle y en el contexto de su naturaleza real. Un modelo matemático
% es la simplificación de ese modelo real para llegar a una solución más
% simple. Requiere de uno o varias hipotesis para reducir la complejidad.
% Por ejemplo, para esta curso se utilizan las hipotesis de que los
% materiales son isotropos y homegeneos.
%% EJ 3
clear; close all; clc
%% Datos
b = 50;
h = 100;
I = b*h^3/12;
E = 70000;
L = 1000;
q = 1;
M = 1000000/2;
%% Nodos
cnod = [0 0
        L 0
        3*L 0];
nnod = 3;
long = [L 2*L];
elem = [1 2
        2 3];
nelem = 2;
%% Dof
dof = [1 2
       3 4
       5 6];
%% Matriz
Kglobal = zeros(2*3);
for e = 1:nelem
    Y4 = 2*E*I/long(e);
    Y3 = Y4*2;
    Y2 = Y4*3/long(e);
    Y1 = Y2*2/long(e);
    Kviga=[Y1 Y2 -Y1 Y2
           Y2 Y3 -Y2 Y4
          -Y1 -Y2 Y1 -Y2
           Y2 Y4 -Y2 Y3];
    dofr = [dof(elem(e,1),:) dof(elem(e,2),:)];
    Kglobal(dofr,dofr) = Kglobal(dofr,dofr)+Kviga;
end
%% BC
fijo = [1 1 1 0 1 0];
libre = ~fijo;
%% Cargas
P = zeros(nnod*2,1);
for i = 1:nelem
    m = i*2-1;
    n = i*2;
    P(m) = P(m)-q*long(i)/2;
    P(m+2) = P(m+2)-q*long(i)/2;
    P(n) = P(n)-q*long(i)^2/12;
    P(n+2) = P(n+2)+q*long(i)^2/12;
end
P(6) = P(6)-M;
%% Reducción
Kglobalred = Kglobal(libre,libre);
Pred = P(libre);
%% Solver
Dred = Kglobalred\Pred;
D = [0 0 0 Dred(1) 0 Dred(2)]';
Dcompleta = [0 0
             0 Dred(1)
             0 Dred(2)
             0 Dred(1)
             0 0]
%% Capor de desplazamientos
X1 = linspace(0,L,20);
X2 = linspace(0,2*L,40);
X = [X1 X2];
e=1;
N1=@(x) 1-3*(x/long(e)).^2+2*(x/long(e)).^3;
N2=@(x) x-2*x.^2/long(e)+x.^3/long(e)^2;
N3=@(x) 3*(x/long(e)).^2-2*(x/long(e)).^3;
N4=@(x) -x.^2/long(e)+x.^3/long(e)^2;
Ne1 = [N1(X(1:19)); N2(X(1:19)); N3(X(1:19)); N4(X(1:19))];
e=2;
N1=@(x) 1-3*(x/long(e)).^2+2*(x/long(e)).^3;
N2=@(x) x-2*x.^2/long(e)+x.^3/long(e)^2;
N3=@(x) 3*(x/long(e)).^2-2*(x/long(e)).^3;
N4=@(x) -x.^2/long(e)+x.^3/long(e)^2;
Ne2 = [N1(X(21:60)); N2(X(21:60)); N3(X(21:60)); N4(X(21:60))];
Desp = [D(1:4)'*Ne1 D(3:6)'*Ne2]';