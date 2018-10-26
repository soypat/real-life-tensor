%Q8, tau_xy
TestType ='Tau_xy';% 'Sigma_y' ; % 'Sigma_y' ; % 'Sigma_x' ; %
eleType = 'Q8';% 
escala=1; %Escala desplazamientos
%% Discretizacion
nodos = [2 0;0 2;0 1;1 0;2*sind(45) 2*sind(45);0 1.5;sind(45) sind(45);1.5 0];
elementos = 1:8;

%% Propiedades del Material {plain stress}
E=207e3;
nu=0.3;
% C = (E/(1-nu^2))*[1 nu 0;
%                   nu 1 0;
%                   0 0 (1-nu)/2];
%%Plane Strain
C = (E/((1+nu)*(1-2*nu)))*[1-nu    nu      0;
                           nu  1-nu      0;
                            0    0  0.5-nu];

%% Definiciones    (N=Número)
Ndofpornod = 2;                    % grados de libertad por nodo
Nelem = size(elementos,1);         % elementos
Nnod = size(nodos,1);           % nodos
Nnodporelem = size(elementos,2);     % nodos por elemento
doftot = Ndofpornod*Nnod;         % grados de libertad
Ndims = size(nodos,2);          % dimensiones del problema
%% Condiciones de borde
bc = false(Nnod,Ndofpornod);       % Matriz de condiciones de borde
bc([4 8 1],2) = true;
bc([3 6 2],1) = true;
%%
% figure(1)
% Meshplot(elementos,nodos,bc,'k',1)

%% Cargas
R = zeros(Nnod,Ndofpornod);
R = zeros(doftot,1);% Vector de cargas

rsInt = 3*ones(1,2); %2 puntos de gauss. Recibe cuantos 
[wpg, gp] = gauss1D(3);
upg=[-sqrt(.6) -1
     0 -1
     sqrt(.6) -1];
npg=3;
%% Carga a Presión
nodosSuperficie=[1 5 2];
q=10;
for iele = 1:Nelem
eleDofs = node2dof(elementos(iele,:),Ndofpornod);
meindof=reshape(eleDofs,2,[])';
for ipg = 1:npg %ultimo punto de gauss
    % Punto de Gauss
    ksi = upg(ipg,1);
    eta = upg(ipg,2);
    N = shapefuns([ksi eta],eleType);

    jac = dN*nodesEle;

    dNxy = jac\dN;          % dNxy = inv(jac)*dN

    for snod=nodosSuperficie %OJO CON ESTA PARTE. Itero SOBRE la matriz. es cosa de pros.
        %snod=nodosSuperficie(i);
        R(meindof(snod,1)) = R(meindof(snod,1))-N(snod)*q*jac(1,2)*wpg(ipg);  %fza en x sobre nodo 2 Q8
        R(meindof(snod,2)) = R(meindof(snod,2))+N(snod)*q*Jac(1,1)*wpg(ipg);  %fza en y sobre nodo 2 Q8
    end
end  
end

%% Puntos de Gauss Matriz rigidez
rsInt = 3*ones(1,2); %2 puntos de gauss. Recibe cuantos 
[wpg, upg, npg] = gauss(rsInt);

%% Matriz de rigidez
K = zeros(doftot);
A = 0; %
jmin = 1E10;%guardamos j min for our information
for iele = 1:Nelem
    Ke = zeros(Ndofpornod*Nnodporelem);
    nodesEle = nodos(elementos(iele,:),:);
    for ipg = 1:npg %ultimo punto de gauss
        % Punto de Gauss
        ksi = upg(ipg,1);
        eta = upg(ipg,2);
        
        % Funciones de forma respecto de ksi, eta
        N = shapefuns([ksi eta],eleType);
        %shapefuns([-1 1],eleType); nos da [0 0 0 1] porque estamos parado
        %en el lado superior izquierdo del elemento
        % Derivadas de las funciones de forma respecto de ksi, eta
        dN = shapefunsder([ksi eta],eleType);
        %sum(dN,2) siempre da uno
        % Derivadas de x,y, respecto de ksi, eta
        jac = dN*nodesEle; %Mucho lio para UN punto del elemento
        %dN es el mismo para el elemento, lo que cambia es nodesEle para
        %obtener el jabociano
        % Derivadas de las funciones de forma respecto de x,y.
        dNxy = jac\dN;          % dNxy = inv(jac)*dN
        
        B = zeros(size(C,2),Ndofpornod*Nnodporelem);
        B(1,1:2:Ndofpornod*Nnodporelem-1) = dNxy(1,:);
        B(2,2:2:Ndofpornod*Nnodporelem) = dNxy(2,:);
        B(3,1:2:Ndofpornod*Nnodporelem-1) = dNxy(2,:);
        B(3,2:2:Ndofpornod*Nnodporelem) = dNxy(1,:);
        
        Djac = det(jac);
        Ke = Ke + B'*C*B*wpg(ipg)*Djac;
        A = A + wpg(ipg)*Djac;
        if Djac < jmin
            jmin = Djac;
        end
    end
    eleDofs = node2dof(elementos(iele,:),Ndofpornod);
    K(eleDofs,eleDofs) = K(eleDofs,eleDofs) + Ke;
end

% Reduccion Matriz
isFixed = reshape(bc',[],1);
isFree = ~isFixed;

Rr = reshape(R',[],1);

% Solver
Dr = K(isFree,isFree)\Rr(isFree);

% Reconstrucción
D = zeros(doftot,1);
D(isFree) = D(isFree) + Dr;

% Reacciones
Rv = K(isFixed,isFree)*D(isFree);
reacciones = nan(doftot,1);
reacciones(isFixed) = Rv;
reacciones = (reshape(reacciones,Ndofpornod,[]))';

%% Recuperación de tensiones en los nodos
uNod = [-1 -1
         1 -1
         1  1
        -1  1
         0 -1
         1  0
         0  1
        -1  0];%Q8 uNod

stress = zeros(Nnodporelem,3,Nelem);
for iele = 1:Nelem
    nodesEle = nodos(elementos(iele,:),:);
    for inode = 1:Nnodporelem
        % Punto de Gauss
        ksi = uNod(inode,1);
        eta = uNod(inode,2);
        % Derivadas de las funciones de forma respecto de ksi, eta
        dN  = shapefunsder([ksi eta],eleType);
        % Derivadas de x,y, respecto d % '2'; % e ksi, eta
        jac  = dN*nodesEle;
        % Derivadas de las funciones de forma respecto de x,y.
        dNxy = jac\dN;          % dNxy = inv(jac)*dN
        
        B = zeros(size(C,2),Ndofpornod*Nnodporelem);
        B(1,1:2:Ndofpornod*Nnodporelem-1) = dNxy(1,:);
        B(2,2:2:Ndofpornod*Nnodporelem) = dNxy(2,:);
        B(3,1:2:Ndofpornod*Nnodporelem-1) = dNxy(2,:);
        B(3,2:2:Ndofpornod*Nnodporelem) = dNxy(1,:);
        
        eleDofs = node2dof(elementos(iele,:),Ndofpornod);
        stress(inode,:,iele) = C*B*D(eleDofs);% Aca esta resuelto el problem
        %C*B es sigma, por los desplazaments
    end
end
    
%% Configuración deformada
nodePosition = nodos + escala*(reshape(D,Ndofpornod,[]))';
figure(1)
Meshplot(elementos,nodePosition,bc,'r',0)
%% Gráfico
figure(2)
bandplot(elementos,nodePosition,stress(:,S,:)',[],'k');