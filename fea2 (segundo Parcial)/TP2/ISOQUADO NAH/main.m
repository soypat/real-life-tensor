%% Comienza la epica
aux=load('TPele.txt');
elementos=aux(:,2:9);%porque es Q4
aux=load('TPnod.txt');
nodos=aux(:,2:3);
nodenumbering=aux(:,1);
nodxs=-1*ones(max(nodenumbering),2); %Creo una nueva matriz que tiene los nodos en la fila correspondiente a su numero asignado
for i=1:size(nodos,1)
    nodxs(nodenumbering(i),[1 2])=nodos(i,:); %Tengo que hacer esto porque Mati o quien sea programo esta porqueria para que tome los nodos segun su linea. Parece a proposito este engendro de paradigma de programacion, seguro lo es.
end
fakeNodeNumber=max(nodenumbering);
fNN=fakeNodeNumber;
nodos=nodxs; %Finalmente asigno la matriz nodos
% Seguro podria hacer el programa de cero, pero para que? No quiero
% aprender.
E=1; 
nu=0.3;
C = (E/(1-nu^2))*[1 nu 0;
                  nu 1 0;
                  0 0 (1-nu)/2]; %Plain Stress
          % tal cual como me lo dieron a mi:
nDofNod = 2;                    % grados de libertad por nodo
nel = size(elementos,1);         % elementos
nNod = size(nodos,1);           % nodos
nNodEle = size(elementos,2);     % nodos por elemento
nDofTot = nDofNod*nNod;         % grados de libertad
nDims = size(nodos,2);          % dimensiones del problema
%% CB
bc = false(nNod,nDofNod);       % Matriz de condiciones de borde

fixity
%% PLOT
figure(1)
myMeshplot(elementos,nodos,bc,'k',1,0) %Eligo si quiero con numeración con el ultimo parametro (1/0)

%% Cargas
R = zeros(nNod,nDofNod);
% %Fuerzas de volumen (Fuerza centrifuga!!!!)
% fv=reshape([zeros(1,nNod);-9.81*2000*ones(1,nNod)],nDofTot,1);
% Fv=int*fv; Fv=reshape(Fv,nDofNod,nNod)';
% R=R+Fv;

%Carga Lineal
% Fuerzas en superficies
q1 =@(x) .01*(150 - (210-x))/150;
% q2 = @(x) -550000*(1-x/52.5);
wallnodes=load('loadnod.txt');
wnodespos=nodos(wallnodes,2);
Lwall=diff(wnodespos);
Nwall=size(wallnodes,1);
Nq8=[4 2 -1;2 16 2;-1 2 4];
Fv=@(qv,L) L/15*Nq8*qv;

for iele=1:Nwall
    
    if iele<=(Nwall-1)/2 %para que no salga de indice
        index=[iele*2-1 iele*2 iele*2+1];
        pos=wnodespos(index);
        qv=[q1(pos(1)); q1(pos(2)) ;q1(pos(3))];
        rvl=(Lwall(iele*2,1)/15)*Nq8*qv;
        R(wallnodes(iele*2-1,1),1)=R(wallnodes(iele*2-1,1),1)+rvl(1);
        R(wallnodes(iele*2,1),1)=R(wallnodes(iele*2,1),1)+rvl(2);
        R(wallnodes(iele*2+1,1),1)=R(wallnodes(iele*2+1,1),1)+rvl(3);
        
    end
    %No esta definido para otro tipo de elemento que no sea Q8
end
%% Puntos de Gauss
rsInt = 3*ones(1,2);
[wpg, upg, npg] = gauss(rsInt);

%% Matriz de rigidez
K = zeros(nDofTot);
A = 0;
jmin = 1E10;
int = zeros(nDofTot);

for iele = 1:nel
    
    Ke = zeros(nDofNod*nNodEle);
    nodesEle = nodos(elementos(iele,:),:);
    eleDofs = node2dof(elementos(iele,:),nDofNod);
    eleDofs = reshape(eleDofs',[],1);
    
    for ipg = 1:npg
        % Punto de Gauss
        ksi = upg(ipg,1);
        eta = upg(ipg,2);
        
        % Funciones de forma respecto de ksi, eta
        N = shapefuns([ksi eta],eleType);
        
%         Nm(1,1:2:7)=N;
%         Nm(2,2:2:8)=N;
        
        % Derivadas de las funciones de forma respecto de ksi, eta
        dN = shapefunsder([ksi eta],eleType);
        % Derivadas de x,y, respecto de ksi, eta
        jac = dN*nodesEle;
        % Derivadas de las funciones de forma respecto de x,y.
        dNxy = jac\dN;          % dNxy = inv(jac)*dN
        
        B = zeros(size(C,2),nDofNod*nNodEle);
        B(1,1:2:nDofNod*nNodEle-1) = dNxy(1,:);
        B(2,2:2:nDofNod*nNodEle) = dNxy(2,:);
        B(3,1:2:nDofNod*nNodEle-1) = dNxy(2,:);
        B(3,2:2:nDofNod*nNodEle) = dNxy(1,:);
        
        Djac = det(jac);
        Ke = Ke + B'*C*B*wpg(ipg)*Djac;
        A = A + wpg(ipg)*Djac;
        %hago la integral me va a servir despues para calcular las fuerzas
%         int(eleDofs,eleDofs)=int(eleDofs,eleDofs)+Nm'*Nm*det(jac);
        if Djac < jmin
            jmin = Djac;
        end
    end
    
    K(eleDofs,eleDofs) = K(eleDofs,eleDofs) + Ke;
end
