%% Elemento exotico Q5 rectangular 2D. (Q4 con nodo extra en 0,0)
% clear; close all; clc
% clear
%% Obtención de funciones de forma
syms x y real %Q5
X = [1 x y x^2 x*y];
A = [1 -1 -1  1  1
     1  1 -1  1 -1
     1  1  1  1  1
     1 -1  1  1 -1
     1  0  0  0 0];  
N = X/A;
dN = [diff(N,x); diff(N,y)];
sum(N);  % == 1
sum(dN,2); % == 0

Ncook(1,1:2:2*length(N))=N;
Ncook(2,2:2:2*length(N))=N; %Tiene la forma de las funciones de forma encontradas en el cook pg 206, ecuacion (6.2-2). Despues veo si me sirven
dNcook(1,1:2:2*length(N))=diff(N,x);
dNcook(2,2:2:2*length(N))=diff(N,y);

%% Datos
E = 200E3; %[MPa]
nu = 0.3;
lam = E*nu/(1+nu)/(1-2*nu);
mu = E/2/(1+nu);
t = 1;

%% Nodos y Elementos
nod = [-1 -1
        1 -1
        1  1
       -1  1
        0  0]*1E3;
nnod = size(nod,1);
elem = [1 2 3 4 5];
nelem = size(elem,1);

%% DOF
dofpornodo = 2;
nodporelem = 5;
doftot = dofpornodo*nnod;
dof = reshape((1:doftot)',dofpornodo,nnod)';

%% Constitutiva (PlaneStrain)
Cstrain = E/((1 + nu)*(1 - 2*nu))*[ 1 - nu      nu            0.0
                                      nu       1 - nu         0.0
                                      0.0        0.0     (1 - 2*nu)/2 ];

%% Gauss para regla de 2x2 (numeración de nodos como figura 6.3-3 pág 212)
a   = 1/sqrt(3);
% Ubicaciones puntos de Gauss
upg = [ -a  -a
        -a   a
         a  -a
         a   a ];
% Número de puntos de Gauss
npg = size(upg,1);
wpg = ones(npg,1); %Weight vale 1 para orden n=2 tabla pg 210.

%% Gauss para regla 3x3
a3 = sqrt(.6);
upg3 = [-a3 -a3; %Para orden n=3 (no lo usamos en este caso)
        -a3  0;
        -a3  a3;
         0  -a3;
         0   0; % Nodo 5 está en el medio, va tener un weight diferente
         0   a3;
         a3  -a3;
         a3   0;
         a3   a3]; %Numerados según figura 6.3-3

wpg3    = 5/9*ones(9,1); %sale de la misma tabla (6.3-1)
wpg3(5) = 8/9;
npg3=9;

%% Matriz de rigidez
Kglobal = zeros(doftot);
for e = 1:nelem
    Ke = zeros(dofpornodo*nodporelem);
    nodelem = nod(elem(e,:),:);
    signo=NaN;
    for ipg = 1:npg
        % Punto de Gauss
        x = upg(ipg,1);
        y = upg(ipg,2);  

        dNs=subs(dN);
        % Derivadas de x,y, respecto de ksi, eta
        
        J = dNs*nodelem;
        if ~isnan(signo) && sign(det(J))~=signo || signo==0 %verifica el bienestar de los elementos
            warning('\nJacobian sign change or singularity detected element %0.0f\n',e)
        end
        signo=sign(det(J));
        
        % Derivadas de las funciones de forma respecto de x,y.
        dNxy = J\dNs; %Primera Forma de calcular B
        dNxycook=J\subs(dNcook); %segunda Forma (cook)

        B = zeros(size(Cstrain,2),size(Ke,1));%Forma 1
        B(1,1:2:9)  = dNxy(1,:);
        B(2,2:2:10) = dNxy(2,:);
        B(3,1:2:9)  = dNxy(2,:);
        B(3,2:2:10) = dNxy(1,:); 
        
        B2=zeros(size(Cstrain,2),size(Ke,1)); %Segunda Forma
        B2(1:2,:)=dNxycook;
        B2(3,2:end)=dNxycook(1,1:end-1);
        B2(3,1:end-1)=B2(3,1:end-1)+dNxycook(2,2:end);

        Ke = Ke + B'*Cstrain*B*wpg(ipg)*det(J);
    end
    dofs = dof(elem(e,:),:);
    dofs = reshape(dofs',[],1);
    Kglobal(dofs,dofs) = Kglobal(dofs,dofs) + Ke;  
end

%% BC
fijo = zeros(doftot/2,2);
fijo([1 2],:) = [1 1; 0 1];
fijo = logical(reshape(fijo',[],1));
libre = ~fijo;

%% Cargas
P = zeros(doftot/2,2);
Q = 0.002; %[N/mm]
f =@(x) Q*x;
a   = 1/sqrt(3);
% a = sqrt(0.6);
upg = [-a a];    
npg = size(upg,2);
wpg = ones(npg,1);
% wpg = [5 8 5]/9;
for e = 1:nelem
    nodelem = nod(elem(e,:),:);
    for ipg = 1:npg
        % Punto de Gauss
        x = upg(ipg);
        y = 1;
        Ns=subs(N);

        dNs=subs(dN);
        sig = f(x);
        J = dNs*nodelem;
        P(elem(e,4),:) = P(elem(e,4),:) + Ns(4)*wpg(ipg)*t*[0 sig*J(1,1)];
        P(elem(e,3),:) = P(elem(e,3),:) + Ns(3)*wpg(ipg)*t*[0 sig*J(1,1)];
    end
end
P = reshape(P',[],1);

%% Solver
Dred = Kglobal(libre,libre)\P(libre);
D = zeros(doftot,1);
D(libre) = Dred;
cnodfinal = nod+reshape(D,dofpornodo,nnod)'*10000000;
figure
hold on; grid on; axis equal;
plot(nod(:,1),nod(:,2),'.')
plot(cnodfinal(:,1),cnodfinal(:,2),'.')

%% Gauss para regla de 2x2 (numeración de nodos como figura 6.3-3 pág 212)
a   = 1/sqrt(3);
% Ubicaciones puntos de Gauss
upg = [ -a  -a
        -a   a
         a  -a
         a   a ];    
% Número de puntos de Gauss
npg = size(upg,1);
wpg = ones(npg,1);

%% Tensiones
stress = zeros(nelem,nodporelem-1,3);
sigvm = zeros(nelem,nodporelem-1);
for e = 1:nelem
    nodelem = nod(elem(e,1:4),:);
    for in = 1:npg
        x = upg(in,1);
        y = upg(in,2);  
%         dN = 1/4*[-(1-y)   1-y    1+y  -(1+y)
%                   -(1-x) -(1+x)   1+x    1-x ];
        dN = [dNx(1:4); dNy(1:4)];
        dNs = subs(dN);
        J = dNs*nodelem;
        
        dNxy = J\dNs;
        B = zeros(size(Cstrain,2),dofpornodo*npg);
        B(1,1:2:7) = dNxy(1,:);
        B(2,2:2:8) = dNxy(2,:);
        B(3,1:2:7) = dNxy(2,:);
        B(3,2:2:8) = dNxy(1,:);
        dofs = reshape(dof(elem(e,1:4),:)',[],1);
        stress(e,in,:) = Cstrain*B*D(dofs);
            end
end

sigvm=sqrt(stress(:,:,1).^2+stress(:,:,2).^2-stress(:,:,1).*stress(:,:,2)+3*stress(:,:,3).^2)
Po=P
stresso=stress
Do=D
sigo=sigvm