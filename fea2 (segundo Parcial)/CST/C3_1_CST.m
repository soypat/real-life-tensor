%% C3.1 CST
close all; clear; clc
%% Datos
E = 200e3; % [MPa]
nu = 0.3;
alfa = 12e-6; % [1/C]
t = 1; % [mm]

graph=false;
%% Shape functions CST
syms x y b h real
X = [1 x y];
A = [1 0 0
     1 1 0
     1 0 1];
N = X/A;
dNx = diff(N,x);
dNy = diff(N,y);
dN = [dNx; dNy];
disp(N);
disp(dN);
%Generados:
N =@(x,y) [ 1 - y - x, x, y];
dN = [ -1, 1, 0; -1, 0, 1]; %Derivadas constantes
%% Nodos y elementos
nod = [0 0
       10 0
       30 0
       0 10
       30 15
       0 20
       15 20
       30 20];
nnod = size(nod,1);
elem = [1 2 4
        2 3 5
        2 7 4
        2 5 7
        4 7 6
        5 8 7];
nelem = size(elem,1);
nnodelem = size(elem,2);
% Meshplot(elem,nod,'b')
% hold on
%% DOF
ndofnod = 2;
doftot = ndofnod*nnod;
dof = reshape(1:doftot,ndofnod,nnod)';

%% BC
bc = false(nnod,ndofnod);
bc(1,:) = true;
bc([4 6],1) = true;

%% Gauss grado 2 tringular
a   = 1/2;
upg = [a  0
       0  a
       a  a];    
npg = size(upg,1);
wpg = ones(npg,1)/3;

%% Constitutiva
C = E/(1 - nu^2)*[ 1.0      nu        0.0
                   nu       1.0       0.0
                   0.0      0.0    (1 - nu)/2 ];

%% Matriz de rigidez
K = zeros(doftot);
A=0;
for e = 1:nelem
    Ke = zeros(ndofnod*nnodelem);
    nodelem = nod(elem(e,:),:);
    
    for ipg=1:npg
    jac = dN*nodelem;
    dNxy = jac\dN;
    
    B = zeros(size(C,2),ndofnod*nnodelem);
    B(1,1:2:5) = dNxy(1,:);
    B(2,2:2:6) = dNxy(2,:);
    B(3,1:2:5) = dNxy(2,:);
    B(3,2:2:6) = dNxy(1,:); 
    A=wpg(ipg)*det(jac)/2+A;
    Ke = Ke + B'*C*B*det(jac)*t;
    
    end
    eleDofs = dof(elem(e,:),:);
    eleDofs = reshape(eleDofs',[],1);
    K(eleDofs,eleDofs) = K(eleDofs,eleDofs) + Ke;
end

%% Cargas
P = 100000; % [N/m^2];
R = zeros(nnod,ndofnod);
R([3 5],1) = R([3 5],1) + P*15/2;
R([5 8],1) = R([3 8],1) + P*5/2;

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
        0 1];
for e = 1:nelem
    nodelem = nod(elem(e,:),:);
    for in = 1:nnodelem
        jac = dN*nodelem;                      
        % Derivadas de las funciones de forma respecto de x,y.
        dNxy = jac\dN;
        B = zeros(size(C,2),ndofnod*nnodelem);
        B(1,1:2:5) = dNxy(1,:);
        B(2,2:2:6) = dNxy(2,:);
        B(3,1:2:5) = dNxy(2,:);
        B(3,2:2:6) = dNxy(1,:);
        dofs = reshape(dof(elem(e,:),:)',[],1);
        stress(e,in,:) = C*B*D(dofs);
    end
end
%% Configuracion deformada
D = (reshape(D,ndofnod,[]))';
nodePosition = nod + D(:,1:2);

%Gráfico
if graph==true
bandplot(elem,nodePosition,stress(:,:,1),[],'k');
Meshplot(elem,nod,bc,'k',1)

end





























