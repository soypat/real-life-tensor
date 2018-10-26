clear
% clc
close all
eleType='Q4';
escala=1;
aux=load('TPele.txt');
elementos=aux(:,2:5);%porque es Q8
aux=load('TPnod.txt');
nodos=aux(:,2:3); %mm

numeracion=aux(:,1);
[nodos, elementos]=nodekill(nodos,numeracion,elementos);
% nodxs=-1*ones(max(nodenumbering),2); %Creo una nueva matriz que tiene los nodos en la fila correspondiente a su numero asignado
% for i=1:size(nodos,1)
%     nodxs(nodenumbering(i),[1 2])=nodos(i,:); %Tengo que hacer esto porque Mati o quien sea programo esta porqueria para que tome los nodos segun su linea. Parece a proposito este engendro de paradigma de programacion, seguro lo es.
% end
% fakeNodeNumber=max(nodenumbering);
% fNN=fakeNodeNumber;
% nodos=nodxs; %Finalmente asigno la matriz nodos
% Discretizacion
% aux=load('TPele.txt');
% elementos=aux(:,2:9);%porque es Q4
% aux=load('TPnod.txt');
% nodos=aux(:,2:3);

%% Propiedades del Material
k=25;
C=[k 0;0 k];
%% Definiciones
%reescribimos la matriz elementos y nodos, sino lo cargo de adina
          
nDofNodT = 1;                    % grados de libertad por nodo
nelT = size(elementos,1);         % elementos
nNodT = size(nodos,1);           % nodos
nNodEleT = size(elementos,2);     % nodos por elemento
nDofTotT = nDofNodT*nNodT;         % grados de libertad
nDimsT = size(nodos,2);          % dimensiones del problema

%% Condiciones de borde
bc = false(nDofTotT,1);       % Matriz de condiciones de borde
% bc(2:5,1:2) = true;
%%
% figure(1)
% myMeshplot(elementos,nodos,bc,'k',1,1)


%% Puntos de Gauss
rsInt = 3*ones(1,2);
[wpg, upg, npg] = gauss(rsInt);

%% Matriz de rigidez
K = zeros(nDofTotT);
A = 0;
jmin = 1E10;
for iele = 1:nelT
    
    Ke = zeros(nNodEleT);
    nodesEle = nodos(elementos(iele,:),:);
    eleDofs = node2dof(elementos(iele,:),1);
    eleDofs = reshape(eleDofs',[],1);
    
    for ipg = 1:npg
        % Punto de Gauss
        ksi = upg(ipg,1);
        eta = upg(ipg,2);
        
        % Funciones de forma respecto de ksi, eta
        N = shapefuns([ksi eta],eleType);
                
        % Derivadas de las funciones de forma respecto de ksi, eta
        dN = shapefunsder([ksi eta],eleType);
        % Derivadas de x,y, respecto de ksi, eta
        jac = dN*nodesEle;
        % Derivadas de las funciones de forma respecto de x,y.
        dNxy = jac\dN;          % dNxy = inv(jac)*dN
        
        B = zeros(2,4);
        B(1,:) = dNxy(1,:);
        B(2,:) = dNxy(2,:);
        Djac = det(jac);
        Ke = Ke + B'*C*B*wpg(ipg)*Djac;
        A = A + wpg(ipg)*Djac;
        %hago la integral me va a servir despues para calcular las fuerzas
        if Djac < jmin
            jmin = Djac;
        end
    end

    K(eleDofs,eleDofs) = K(eleDofs,eleDofs) + Ke;
end
%% Cargas
Co=[];
Tc=[];
for i=1:nNodT
    xy=nodos(i,:);
    if (xy(1)==40 || xy(1)==70) && xy(2)>=90
        if xy(1)==40
            Tc=[Tc; 0];
        else
            Tc=[Tc; 120];
        end
        Co=[Co i];
    end
end

X=zeros(1,nNodT-length(Co));
aux=1;
for inod=1:size(nodos,1)
    if isempty(find(Co==inod,1))
        X(aux)=inod;
        aux=aux+1;
    end
end
%% Coso donde haces xx cc xc
% X=[2,3,4,7,8,9];
Kxx=K(X,X);
Kxc=K(X,Co);
Q=Kxc*Tc;
Tx=Kxx\-Q;
T=zeros(15,1);
T(X)=Tx;
T(Co)=Tc;

for i=1:nelT
    temp(i,:)= T(elementos(i,:));
%     flujo(i,:)= q(elementos(i,:));
end

figure(2)
bandplot(elementos,nodos,temp)



