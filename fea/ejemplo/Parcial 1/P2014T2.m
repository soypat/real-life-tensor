%% 2014T2
%% EJ 1
%a. La matriz es simétrica ya que los existe una relación lineal entre los
%desplazamientos y las cargas.
%b. La diagonal debe ser positiva ya que cada elemnto de la misma
%respresenta la relación directa entre un grado de libertad y la respectiva
%fuerza que lo afecta. Es fisicamente imposible que una fuerza desplace un
%grado de libertad en sentido opuesto al suyo.
%% EJ 2
% Un modelo matemático intenta representar la un problema físico con
% simplificaciones de la realidad. Aún así, muchas veces las ecuaciones
% resultantes son muy dificiles de resolver. Un modelo numérico es una
% simplificación lineal del problema matemático que, en el mejor de los
% casos, convergerá a la solución excacta del modelo matemático.
%% EJ 3
clear; close all; clc
%% Datos
b = 100;
h = 50;
I = b*h^3/12;
E = 200000;
q = 2;
M = 2E6;
L = 1000;
%% Nodos y elementos
% Se utiliza simetría en la mitad del sistema
cnod = [0 0
        L 0
        L*2.5 0
        L*4 0];
nnod = 4;
elem = [1 2
        2 3
        3 4];
nelem = 3;
%% Dof
dofpornodo = 2;
dof = reshape((1:1:dofpornodo*nnod)',dofpornodo,nnod)';
%% Matriz global
Kglobal = zeros(dofpornodo*nnod);
memlong = zeros(1,nelem);
for e = 1:nelem
    v = cnod(elem(e,2),:) - cnod(elem(e,1),:);
    long = norm(v);
    memlong(e) = long;
    Y4 = 2*E*I/long;
    Y3 = Y4*2;
    Y2 = Y4*3/long;
    Y1 = Y2*2/long;
    Kviga=[Y1 Y2 -Y1 Y2
           Y2 Y3 -Y2 Y4
          -Y1 -Y2 Y1 -Y2
           Y2 Y4 -Y2 Y3];
    dofr = [dof(elem(e,1),:) dof(elem(e,2),:)];
    Kglobal(dofr,dofr) = Kglobal(dofr,dofr)+Kviga;
end
%% BC
fijo = 1&[1 1 1 0 0 0 0 1];
libre = ~fijo;
%% Cargas
Pmatriz = [-q*memlong(1)/2                  -q*memlong(1)^2/12
           -q*memlong(1)/2+q*memlong(2)/2    q*memlong(1)^2/12+q*memlong(2)^2/12
           q*memlong(2)/2+q*memlong(3)/2    -q*memlong(2)^2/12+q*memlong(3)^2/12-M
           q*memlong(3)/2                   -q*memlong(3)^2/12];
P = reshape(Pmatriz',[],1);
%% Reducción
Kglobalred = Kglobal(libre,libre);
Pred = P(libre);
% Solver
Dred = Kglobalred\Pred;
D = zeros(dofpornodo*nnod,1);
D(libre) = Dred;
Dtotal = [D(1:4)' 0 -D(4) 0 0]';
disp('Desplazamientos de grados de libertad')
Dtotal = reshape(Dtotal,2,4)'
%% Desplazamientos
Dpornodo = reshape(D,dofpornodo,nnod)';
n = 10;
Desp = zeros(n,nelem);
for e = 1:nelem
    X = linspace(0,memlong(e),n);
    N1=@(x) 1-3*(x/memlong(e)).^2+2*(x/memlong(e)).^3;
    N2=@(x) x-2*x.^2/memlong(e)+x.^3/memlong(e)^2;
    N3=@(x) 3*(x/memlong(e)).^2-2*(x/memlong(e)).^3;
    N4=@(x) -x.^2/memlong(e) +x.^3/memlong(e)^2;
    N = [N1(X); N2(X); N3(X); N4(X)];
    Desp(:,e) = ([Dpornodo(elem(e,1),:) Dpornodo(elem(e,2),:)]*N)';
end
Desplazamientos = [0; reshape(Desp(2:end,:),[],1)]
plot(0:1:27,Desplazamientos)
%% Tensiones en A
X = 500;
e = 1;
DA = [Dpornodo(1,:) Dpornodo(2,:)];
    N1=@(x) -6/memlong(e)^2+12*x/memlong(e)^3;
    N2=@(x) -4/memlong(e)+6*x/memlong(e)^2;
    N3=@(x) 6/memlong(e)^2-12*x/memlong(e)^3;
    N4=@(x) -2/memlong(e)+6*x/memlong(e)^2;
    B = [N1(X); N2(X); N3(X); N4(X)];
SigmaFlex = h/2*E*DA*B