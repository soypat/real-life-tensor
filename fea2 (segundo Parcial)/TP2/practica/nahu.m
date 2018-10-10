clear
% clc
close all

TestType = 'Sigma_x' ; % 'Sigma_y' ; % 'Sigma_x' ; %
eleType = 'Q4';% 'Q4', 'Q9'
escala=1; %Escala desplazamientos
% Discretizacion
aux=load('elementosCQ4.txt');
elementos=aux(:,2:5);%porque es Q4
aux=load('nodosCQ4.txt');
nodos=aux(:,2:3);
numeracion=aux(:,1);
[nodos, elementos]=nodekill(nodos,numeracion,elementos);

%% Propiedades del Material
k=25;
C=[k 0;0 k];
%% Definiciones
%reescribimos la matriz elementos y nodos, sino lo cargo de adina
          
nDofNod = 2;                    % grados de libertad por nodo
nel = size(elementos,1);         % elementos
nNod = size(nodos,1);           % nodos
nNodEle = size(elementos,2);     % nodos por elemento
nDofTot = nDofNod*nNod;         % grados de libertad
nDims = size(nodos,2);          % dimensiones del problema

% %% Condiciones de borde
bc = false(15,2);       % Matriz de condiciones de borde
bc(2:5,1:2) = true;
%%
% figure(1)
% Meshplot(elementos,nodos,bc,'k',1)


%% Puntos de Gauss
rsInt = 3*ones(1,2);
[wpg, upg, npg] = gauss(rsInt);

%% Matriz de rigidez
K = zeros(15);
A = 0;
jmin = 1E10;
for iele = 1:nel
    
    Ke = zeros(nNodEle);
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

Co=[1,6,11,12 13,14,15];
X=[2,3,4,7,8,9];
Kxx=K(X,X);
Kxc=K(X,Co);
Tc=zeros(length(Co),1);
Tc([1,2])=100;
R=Kxc*Tc;
Tx=Kxx\-R;
T=zeros(15,1);
T(X)=Tx;
T(Co)=Tc;

for i=1:nel
    temp(i,:)= T(elementos(i,:));
%     flujo(i,:)= q(elementos(i,:));
end

figure(2)
bandplot(elementos,nodos,temp)

%primero quiere el qk
% q=zeros(2*8,1);
% for iele = 1:nel
%     
%     Ke = zeros(nNodEle);
%     nodesEle = nodos(elementos(iele,:),:);
%     eleDofs = node2dof(elementos(iele,:),1);
%     eleDofs = reshape(eleDofs',[],1);
%     
%     for ipg = 1:npg
%         % Punto de Gauss
%         ksi = upg(ipg,1);
%         eta = upg(ipg,2);
%         
%         % Funciones de forma respecto de ksi, eta
%         N = shapefuns([ksi eta],eleType);
%                 
%         % Derivadas de las funciones de forma respecto de ksi, eta
%         dN = shapefunsder([ksi eta],eleType);
%         % Derivadas de x,y, respecto de ksi, eta
%         jac = dN*nodesEle;
%         % Derivadas de las funciones de forma respecto de x,y.
%         dNxy = jac\dN;          % dNxy = inv(jac)*dN
%         
%         B = zeros(2,4);
%         B(1,:) = dNxy(1,:);
%         B(2,:) = dNxy(2,:);
%         Djac = det(jac);
%         q([2*iele-1 2*iele])=25*B*T(eleDofs);
% 
%         %hago la integral me va a servir despues para calcular las fuerzas
%         if Djac < jmin
%             jmin = Djac;
%         end
%     end
% end
% figure(2)
% qx=q(1:2:end);
% qy=q(2:2:end);
% quiver(qx,qy)
% bandplot(elementos,nodos,qx,qy,'k');
