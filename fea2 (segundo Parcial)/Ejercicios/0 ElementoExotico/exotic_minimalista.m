%% Elemento exotico Q5 rectangular 2D. (Q4 con nodo extra en 0,0)
%  clear; close all; clc
% reset(symengine)
%% Obtención de funciones de forma
clear dN N dNaux
syms x y real
X = [1 x y x^2 x*y];
A = [1 -1 -1  1  1
     1  1 -1  1 -1
     1  1  1  1  1
     1 -1  1  1 -1
     1  0  0  0 0];  
shapefuns = X/A;

N(1,1:2:2*length(shapefuns))=shapefuns;
N(2,2:2:2*length(shapefuns))=shapefuns; %Tiene la forma de las funciones de forma encontradas en el cook pg 206, ecuacion (6.2-2). Despues veo si me sirven
dN(1,1:2:2*length(shapefuns))=diff(shapefuns,x);
dN(2,2:2:2*length(shapefuns))=diff(shapefuns,y);
u=1;
% dNcook(1,1:2:2*length(shapefuns))=diff(shapefuns,x);
% dNcook(2,2:2:2*length(N))=diff(N,y);

dNaux=[diff(shapefuns,x);diff(shapefuns,y)]; %Para calcular jacobiano
%% Datos
E = 200E3; %[MPa]
nu = 0.3;
% lam = E*nu/(1+nu)/(1-2*nu);
% mu = E/2/(1+nu);
t = 1;

%% Nodos y Elementos
nod = [-1 -1
        1 -1
        0  0
       -1  1
        1  1]*1E3;
nnod = size(nod,1);
elem = [1 2 5 4 3];
nelem = size(elem,1);

%% DOF
dofpornodo = 2;
nodporelem = 5;
doftot = dofpornodo*nnod;
dof = reshape((1:doftot)',dofpornodo,nnod)';

%% Constitutiva (PlaneStrain)
Cstrain = E/((1 + nu)*(1 - 2*nu))*[ 1 - nu      nu            0.0;nu       1 - nu         0.0;0.0        0.0     (1 - 2*nu)/2 ];
%% Gauss para regla de 2x2 (numeración de nodos como figura 6.3-3 pág 212)
a   = 1/sqrt(3);
upg = [ -a  -a; -a   a;a  -a;a   a ];% Ubicaciones puntos de Gauss
npg = size(upg,1);
wpg = ones(npg,1); %Weight vale 1 para orden n=2 tabla pg 210.
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
        
        J =subs(dNaux)*nodelem;
        
        if ~isnan(signo) && sign(det(J))~=signo || signo==0 %verifica el bienestar de los elementos
            warning('\nJacobian sign change or singularity detected element %0.0f\n',e)
        end
        signo=sign(det(J));

        dNxy=J\subs(dN);
        
        B=zeros(size(Cstrain,2),size(Ke,1)); %Segunda Forma
        B(1:2,:)=dNxy;
        B(3,2:end)=dNxy(1,1:end-1);
        B(3,1:end-1)=B(3,1:end-1)+dNxy(2,2:end);
        
        
        Ke = Ke + B'*Cstrain*B*wpg(ipg)*det(J);
    end
    dofs = dof(elem(e,:),:);
    dofs = reshape(dofs',[],1);
    Kglobal(dofs,dofs) = Kglobal(dofs,dofs) + Ke;
end

%% BC
fixity=false(doftot,1);
fixity([1 2 4])=true;
free=~fixity;

%% Cargas
P = zeros(doftot/2,2);
q = 0.002; %[N/mm]
f =@(x) q*x;
a   = 1/sqrt(3);
% a = sqrt(0.6);
upg = [-a a];    
npg = size(upg,2);
wpg = ones(npg,1);
% wpg = [5 8 5]/9;
for e = 1:nelem
    nodelem = nod(elem(e,:),:);
    for ipg = 1:npg
        x = upg(ipg);
        y = 1; %El truco que no te querían mostrar para aplicar una carga superficial. Te parás sobre la linea!
        Ns=subs(N);

        dNs=subs(dN);
        sig = f(x);
        J = subs(dNaux)*nodelem;
        Nedge=subs(N(:,3*2-1:4*2));
        
        P(elem(e,4),:) = P(elem(e,4),:) + subs(shapefuns(4))*wpg(ipg)*[0 sig*J(1,1)]; %La integral 6.9-5 del Cook
        P(elem(e,3),:) = P(elem(e,3),:) + subs(shapefuns(3))*wpg(ipg)*[0 sig*J(1,1)];
    end
end
P = reshape(P',[],1);

%% Solver
Dred = Kglobal(free,free)\P(free);
D = zeros(doftot,1);
D(free) = Dred;
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
        dN = 1/4*[-(1-y)   1-y    1+y  -(1+y)
                  -(1-x) -(1+x)   1+x    1-x ];
        
        J = dN*nodelem;
        
        dNxy = J\dN;
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
bandplot(nod,elem,sigvm')