%% C3.1 CST
close all; clear; clc
%% Datos
E = 200e3; % [MPa]
nu = 0.3;
alfa = 12e-6; % [1/C]
t = 1; % [mm]

%% Shape functions CST
syms x y real
X = [1 x y x^2 x*y y^2];
A = [1 0 0 0 0 0
     1 1 0 1 0 0
     1 0 1 0 0 1
     1 .5 0 .25 0 0
     1 .5 .5 .25 .25 .25
     1 0 .5 0 0 .25];
N = X/A;

dNx = diff(N,x);
dNy = diff(N,y);
dN = [dNx; dNy];
disp(N);
disp(dN);
%Generados:
N =@(x,y) [ 2*x^2 + 4*x*y - 3*x + 2*y^2 - 3*y + 1, 2*x^2 - x, 2*y^2 - y, 4*x - 4*x*y - 4*x^2, 4*x*y, 4*y - 4*x*y - 4*y^2];
dN =@(x,y) [ 4*x + 4*y - 3, 4*x - 1,       0, 4 - 4*y - 8*x, 4*y,          -4*y
             4*x + 4*y - 3,       0, 4*y - 1,          -4*x, 4*x, 4 - 8*y - 4*x];
%% Nodos y elementos
nod = load('nodLST.txt');
nod = nod(:,2:3);
nnod = size(nod,1);
elem = load('elemLST.txt');
elem = elem(:,2:7);
nelem = size(elem,1);
nnodelem = size(elem,2);
%meshplot(elem,nod,'b')

%% DOF
ndofnod = 2;
doftot = ndofnod*nnod;
dof = reshape(1:doftot,ndofnod,nnod)';

%% BC
bc = false(nnod,ndofnod);
bc(1,:) = true;
bc(nod(:,1)==0,1) = true;

%% Gauss grado 2 tringular
a   = 1/2;
upg = [a  0
       0  a
       a  a];    
npg = size(upg,1);
wpg = [1 1 1]/3;

%% Constitutiva
C = E/(1 - nu^2)*[ 1.0      nu        0.0
                   nu       1.0       0.0
                   0.0      0.0    (1 - nu)/2 ];

%% Matriz de rigidez
K = zeros(doftot);
for e = 1:nelem
    Ke = zeros(ndofnod*nnodelem);
    nodelem = nod(elem(e,:),:);
    for ipg = 1:npg
        r = upg(ipg,1);
        s = upg(ipg,2);
        dNrs = dN(r,s);
        jac = dNrs*nodelem;
        dNxy = jac\dNrs;
        
        B = zeros(size(C,2),ndofnod*nnodelem);
        B(1,1:2:11) = dNxy(1,:);
        B(2,2:2:12) = dNxy(2,:);
        B(3,1:2:11) = dNxy(2,:);
        B(3,2:2:12) = dNxy(1,:); 

        Ke = Ke + B'*C*B*wpg(ipg)*det(jac)*t;
    end
    eleDofs = dof(elem(e,:),:);
    eleDofs = reshape(eleDofs',[],1);
    K(eleDofs,eleDofs) = K(eleDofs,eleDofs) + Ke;
end

%% Cargas
P = 100000; % [N/m^2];
R = zeros(nnod,ndofnod);
R([7 10 8],1) = R([7 10 8],1) + P*15*[1/6 2/3 1/6]';
R([8 20 19],1) = R([8 20 19],1) + P*5*[1/6 2/3 1/6]';

%% Reduccion Matriz
isFixed = reshape(bc',[],1);
isFree = ~isFixed;

Rr = reshape(R',[],1);

% Solver
Dr = K(isFree,isFree)\Rr(isFree);

% Reconstrución
D = zeros(doftot,1);
D(isFree) = D(isFree) + Dr;

%% Tensiones en los nodos
stress = zeros(nelem,nnodelem,3);
unod = [0 0
        1 0
        0 1
        .5 0
        .5 .5
        0 .5];
for e = 1:nelem
    nodelem = nod(elem(e,:),:);
    for in = 1:nnodelem
        dNrs = dN(unod(in,1), unod(in,2));
        jac = dNrs*nodelem;                      
        % Derivadas de las funciones de forma respecto de x,y.
        dNxy = jac\dNrs;
        B = zeros(size(C,2),ndofnod*nnodelem);
        B(1,1:2:11) = dNxy(1,:);
        B(2,2:2:12) = dNxy(2,:);
        B(3,1:2:11) = dNxy(2,:);
        B(3,2:2:12) = dNxy(1,:);
        dofs = reshape(dof(elem(e,:),:)',[],1);
        stress(e,in,:) = C*B*D(dofs);
    end
end
%% Configuracion deformada
D = (reshape(D,ndofnod,[]))';
nodePosition = nod + D(:,1:2);

%Gráfico
bandplot(elem,nodePosition,stress(:,:,1),[],'k');
meshplot(elem,nod,'b')



























