%% Codigos prearmados

%% Funciones de forma
syms x y real
X = [1 x y x^2 x*y];
A = [1 -1 -1  1  1
     1  1 -1  1 -1
     1  1  1  1  1
     1 -1  1  1 -1
     1  0  0  0  0];
N = X/A;
dNx = diff(N,x);
dNy = diff(N,y);
dN = [dNx; dNy];
sum(N);  % == 1
sum(dN); % == 0

%% Nodos y elementos
nod = [];
nnod = size(nod,1);
elem = [];
nelem = size(elem,1);
nnodelem = size(elem,2);
% meshplot(elem,nod,'b')

%% DOF
ndofnod = 2;
doftot = ndofnod*nnod;
dof = reshape(1:doftot,ndofnod,nnod)';


%% Matriz constitutiva
% Plane strain
C = [1-nu  nu  0
     nu 1-nu 0
     0 0 (1-2*nu)/2]*E/(1+nu)/(1-2*nu);
% Plane stress
C = [1 nu 0
     nu 1 0
     0 0 (1-nu)/2]*E/(1-nu^2);

%% Gauss1D 2 puntos
a   = 1/sqrt(3);
upg = [-a a];    
npg = size(upg,2);
wpg = ones(npg,1);

%% Gauss1D 3 puntos
a = sqrt(0.6);
upg = [-a 0 a];    
npg = size(upg,2);
wpg = [5 8 5]/9;
%% Gauss grado 2 tringular
a   = 1/2;
upg = [a  0
       0  a
       a  a];    
npg = size(upg,1);
wpg = ones(npg,1)/3;

%% Gauss para regla de 2x2
a   = 1/sqrt(3);
upg = [ -a  -a
        -a   a
         a  -a
         a   a ];    
npg = size(upg,1);
wpg = ones(npg,1);

%% Gauss para regla de 3x3       
a   = sqrt(0.6);
upg = [ -a  -a
        -a   0
        -a   a
         0  -a    
         0   0
         0   a 
         a  -a
         a   0
         a   a ];
npg = size(upg,1);
wpg = [5/9, 5/9, 5/9, 5/9, 8/9, 5/9, 5/9, 5/9, 5/9];

%% Matriz de rigidez Q4
Kglobal = zeros(doftot);
for e = 1:nelem
    nodelem = nod(elem(e,:),:);
    Ke = zeros(nodporelem*dofpornodo);
    for ipg = 1:npg
        ksi = upg(ipg,1);
        eta = upg(ipg,2);
        
        dNke = shapefunsder([ksi eta],eleT);
        J = dNke*nodelem;
        dNxy = J\dNke;
        
        B = zeros(size(C,2),nodporelem*dofpornodo);
        B(1,1:2:7) = dNxy(1,:);
        B(2,2:2:8) = dNxy(2,:);
        B(3,1:2:7) = dNxy(2,:);
        B(3,2:2:8) = dNxy(1,:);
        
        Ke = Ke + B'*C*B*t*wpg(ipg)*det(J);
    end
    dofs = reshape(dof(elem(e,:),:)',[],1);
    Kglobal(dofs,dofs) = Kglobal(dofs,dofs) + Ke;
end

%% BC
fijo = zeros(nnod,ndofnod);
fijo = logical(reshape(fijo',[],1));
libre = ~fijo;

%% Solver
Dred = Kglobal(libre,libre)\P(libre);

D = zeros(doftot,1);
D(libre) = Dred;

nodfinal = nod+reshape(D,dofpornodo,nnod)'*1000;
figure
hold on; grid on; axis equal;
plot(nod(:,1),nod(:,2),'.')
plot(nodfinal(:,1),nodfinal(:,2),'.')

%% Cargas volumétricas
R = zeros(nNod,nDofNod);
for e = 1:nel
    nodesEle = nodes(elements(e,:),:);
    for ipg = 1:npg
        ksi = upg(ipg,1);
        eta = upg(ipg,2);  
        N = [];
        dN = []; 
        jac = dN*nodesEle;                      
        R(elements(e,:),2) = R(elements(e,:),2) + N'*rho*g*t*wpg(ipg)*det(jac);
    end
end

%% Cargas superficiales
q =@(y) 1; % [N/m^2]
Q = integral(q,0,1);
qcheck = 0;

for iele = 1:nel
    nodesEle = nodes(elements(iele,:),:);
    if nodesEle(4,2)==60
        for ipg = 1:npg
            ksi = 1;
            eta = upg(ipg);  
            N = [];
            dN = [];  
            jac = dN*nodesEle;
            N = [N(2) N(3)];
            Q = [q(nodesEle(2,2)) q(nodesEle(3,2))]';
            R(elements(iele,2:3),1) = R(elements(iele,2:3),1) ...
                + N'*N*t*wpg(ipg)*Q*jac(2,2);
            qcheck = qcheck + sum(N'*N*t*wpg(ipg)*Q*jac(2,2));
        end
    end
end


% Evaluando en puntos de gauss estructurales
    if nodesEle(3,1)==30&nodesEle(3,2)<=52.5
        for ipg = 1:npg
            ksi = 1;
            eta = upg(ipg);  
            N = [];
            dN = [];  
            delt = (nodesEle(3,2)-nodesEle(2,2))/2;
            cennod = (nodesEle(3,2)+nodesEle(2,2))/2;
            jac = dN*nodesEle;                      
            R(elements(iele,2),1) = R(elements(iele,2),1) + ...
                N(2)*t*wpg(ipg)*q2(cennod+delt*eta)*jac(2,2);
            R(elements(iele,3),1) = R(elements(iele,3),1) + ...
                N(3)*t*wpg(ipg)*q2(cennod+delt*eta)*jac(2,2);
            qcheck2 = qcheck2 + ...
                (N(2)+N(3))*t*wpg(ipg)*q2(cennod+delt*eta)*jac(2,2);
        end
    end
%% Tensiones
stress = zeros(nelem,nodporelem,3);
unod = [ -1 -1
          1 -1
          1  1
         -1  1 ];
for e = 1:nelem
    nodelem = nod(elem(e,:),:);
    for in = 1:nodporelem
        ksi = unod(in,1);
        eta = unod(in,2);  
        dNke = [];
        J = dNke*nodelem;                      
        dNxy = J\dNke;
        B = zeros(size(C,2),dofpornodo*nodporelem);
        B(1,1:2:7) = dNxy(1,:);
        B(2,2:2:8) = dNxy(2,:);
        B(3,1:2:7) = dNxy(2,:);
        B(3,2:2:8) = dNxy(1,:);
        dofs = reshape(dof(elem(e,:),:)',[],1);
        stress(e,in,:) = C*B*D(dofs);
    end
end
% En lo puntos de Gauss es igual pero evaluado en upg

bandplot(elem,nodfinal,stress(:,:,1),[],'k')