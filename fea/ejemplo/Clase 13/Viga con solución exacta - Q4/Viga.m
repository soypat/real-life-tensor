close all
clear all
clc

Nmallas=5;

Dmax=zeros(1,Nmallas);
E_ene=zeros(1,Nmallas);
maxvex=zeros(1,Nmallas);
for m=1:Nmallas
    
    switch m
        case 1
            nodes=load('nodos10x2.txt');
            elements=load('elementos10x2.txt');
        case 2
            nodes=load('nodos20x4.txt');
            elements=load('elementos20x4.txt');
        case 3
            nodes=load('nodos30x6.txt');
            elements=load('elementos30x6.txt');
        case 4
            nodes=load('nodos40x8.txt');
            elements=load('elementos40x8.txt');
        case 5
            nodes=load('nodos60x10.txt');
            elements=load('elementos60x10.txt');
    end
    nodes=nodes(:,[2,3]);
    elements=elements(:,2:5);
    
    nDofNod = 2;                    % Nï¿½mero de grados de libertad por nodo
    nNodEle = 4;                    % Nï¿½mero de nodos por elemento
    nel = size(elements,1);         % Nï¿½mero de elementos
    nNod = size(nodes,1);           % Nï¿½mero de nodos
    nDofTot = nDofNod*nNod;         % Nï¿½mero de grados de libertad

    
    % Propiedades del material y geometría
    w = 1;
    c=10;
    L=100;
    E = 1000;
    NU = 0.25;
    I = w*(2*c)^3/12;
    G = E/(2+2*NU);
    P = 80;
        
    % Parámetro de malla
    h(m) = sqrt((L*2*c)/nel);
    
    escala=1;
    Uex=ex_sol(nodes(:,1),nodes(:,2),'fun');
    maxUex(m)=max(max(abs(Uex)));     %me guardo el maximo desplazamiento en y exacto
    figure
    meshplot(elements,nodes+escala*Uex,'b')
    %% Matriz Constitutiva (plane stress)
    
    C = E/(1 - NU^2)*[ 1.0     NU         0.0
                        NU    1.0         0.0
                       0.0    0.0     (1 - NU)/2 ];
    
    C_1 = inv(C);
    
    %% Gauss
    a   = 1/sqrt(3);
    % Ubicaciones puntos de Gauss
    upg = [ -a  -a
        a  -a
        a   a
        -a   a ];
    % Nï¿½mero de puntos de Gauss
    npg = size(upg,1);
    wpg = ones(npg,1);
    
    %% Matriz de rigidez
    K = zeros(nDofTot);
    nodeDofs = reshape(1:nDofTot,nDofNod,nNod)';
    for iele = 1:nel
        Ke = zeros(nDofNod*nNodEle);
        nodesEle = nodes(elements(iele,:),:);
        for ipg = 1:npg
            % Punto de Gauss
            ksi = upg(ipg,1);
            eta = upg(ipg,2);
            % Derivadas de las funciones de forma respecto de ksi, eta
            dN = 1/4*[-(1-eta)   1-eta    1+eta  -(1+eta)
                -(1-ksi) -(1+ksi)   1+ksi    1-ksi ];
            % Derivadas de x,y, respecto de ksi, eta
            jac = dN*nodesEle;
            % Derivadas de las funciones de forma respecto de x,y.
            dNxy = jac\dN;          % dNxy = inv(jac)*dN
            
            B = zeros(size(C,2),nDofNod*nNodEle);
            B(1,1:2:7) = dNxy(1,:);
            B(2,2:2:8) = dNxy(2,:);
            B(3,1:2:7) = dNxy(2,:);
            B(3,2:2:8) = dNxy(1,:);
            
            Ke = Ke + B'*C*B*wpg(ipg)*det(jac)*w;
        end
        eleDofs = nodeDofs(elements(iele,:),:);
        eleDofs = reshape(eleDofs',[],1);
        K(eleDofs,eleDofs) = K(eleDofs,eleDofs) + Ke;
    end
    
    %% Vector de cargas
    R = zeros(nDofTot,1);
    
    % Solo corte en borde izquierdo: La carga es {0;-txy}
    LeftBoundary = find(abs(nodes(:,1))<1e-9);
    
    % Corte y tensión normal en borde derecho: La carga es {-sigmax;txy}
    RightBoundary = find(abs(nodes(:,1)-L)<1e-9);
    
    for iele = 1:nel
        Re = zeros(nNodEle,nDofNod);
        nodesEle = nodes(elements(iele,:),:);
        if sum(ismember(elements(iele,:),LeftBoundary)) > 1 || sum(ismember(elements(iele,:),RightBoundary)) > 1
            for e=1:4
                if ismember(elements(iele,e),LeftBoundary) && ismember(elements(iele,mod(e,nNodEle)+1),LeftBoundary)
                    l = norm(nodesEle(1,:)-nodesEle(e,:));
                    v = (nodesEle(1,:)-nodesEle(e,:))/l;
                    X = nodesEle(e,:);
                    X1 = X;
                    X2 = X + v*l/3;
                    X3 = X + v*2*l/3;
                    X4 = X + v*l;
                    
                    S1 = ex_sol(X1(1),X1(2),'stress');
                    S2 = ex_sol(X2(1),X2(2),'stress');
                    S3 = ex_sol(X3(1),X3(2),'stress');
                    S4 = ex_sol(X4(1),X4(2),'stress');
                    
                    % Para la función de forma en e
                    N11 = 1;
                    N12 = 2/3;
                    N13 = 1/3;
                    N14 = 0;
                    
                    f11y = N11*(-S1(3)*w);
                    f12y = N12*(-S2(3)*w);
                    f13y = N13*(-S3(3)*w);
                    f14y = N14*(-S4(3)*w);
                    
                    % Para la función de forma en e+1
                    N21 = 0;
                    N22 = 1/3;
                    N23 = 2/3;
                    N24 = 1;
                    f21y = N21*(-S1(3)*w);
                    f22y = N22*(-S2(3)*w);
                    f23y = N23*(-S3(3)*w);
                    f24y = N24*(-S4(3)*w);
                    
                    %-%
                    f1x = 0;
                    f2x = 0;
                    
                    f1y = (l/8)*(f11y + 3*f12y + 3*f13y + f14y); % Regla Simpson 3/8
                    f2y = (l/8)*(f21y + 3*f22y + 3*f23y + f24y); % Regla Simpson 3/8
                    
                elseif ismember(elements(iele,e),RightBoundary) && ismember(elements(iele,mod(e,nNodEle)+1),RightBoundary)
                    l = norm(nodesEle(mod(e,nNodEle)+1,:)-nodesEle(e,:));
                    v = (nodesEle(mod(e,nNodEle)+1,:)-nodesEle(e,:))/l;
                    X = nodesEle(e,:);
                    X1 = X;
                    X2 = X + v*l/3;
                    X3 = X + v*2*l/3;
                    X4 = X + v*l;
                    
                    S1 = ex_sol(X1(1),X1(2),'stress');
                    S2 = ex_sol(X2(1),X2(2),'stress');
                    S3 = ex_sol(X3(1),X3(2),'stress');
                    S4 = ex_sol(X4(1),X4(2),'stress');
                    
                    % Para la función de forma en e
                    N11 = 1;
                    N12 = 2/3;
                    N13 = 1/3;
                    N14 = 0;
                    
                    f11x = N11*(-S1(1)*w);
                    f12x = N12*(-S2(1)*w);
                    f13x = N13*(-S3(1)*w);
                    f14x = N14*(-S4(1)*w);
                    
                    f11y = N11*(S1(3)*w);
                    f12y = N12*(S2(3)*w);
                    f13y = N13*(S3(3)*w);
                    f14y = N14*(S4(3)*w);
                    
                    % Para la función de forma en e+1
                    N21 = 0;
                    N22 = 1/3;
                    N23 = 2/3;
                    N24 = 1;
                    
                    f21x = N21*(-S1(1)*w);
                    f22x = N22*(-S2(1)*w);
                    f23x = N23*(-S3(1)*w);
                    f24x = N24*(-S4(1)*w);
                    
                    f21y = N21*(S1(3)*w);
                    f22y = N22*(S2(3)*w);
                    f23y = N23*(S3(3)*w);
                    f24y = N24*(S4(3)*w);
                    
                    %-%
                    f1x = -(l/8)*(f11x + 3*f12x + 3*f13x + f14x); % Regla Simpson 3/8
                    f2x = -(l/8)*(f21x + 3*f22x + 3*f23x + f24x); % Regla Simpson 3/8
                    
                    f1y = (l/8)*(f11y + 3*f12y + 3*f13y + f14y); % Regla Simpson 3/8
                    f2y = (l/8)*(f21y + 3*f22y + 3*f23y + f24y); % Regla Simpson 3/8
                else
                    f1x = 0;
                    f2x = 0;
                    f1y = 0;
                    f2y = 0;
                end
                Re([e,mod(e,nNodEle)+1],1:2) = Re([e,mod(e,nNodEle)+1],1:2)+[f1x,f1y;f2x,f2y];
            end
        end
        eleDofs = nodeDofs(elements(iele,:),:);
        eleDofs = reshape(eleDofs',[],1);
        R(eleDofs) = R(eleDofs) + reshape(Re',[],1);
    end
    
    %%  Matriz de condiciones de borde
    bc = false(nNod,nDofNod);
    a=find(abs(nodes(:,1)-L)<1e-8);
    b1=find(abs(nodes(a,2)-c)<1e-8);
    b2=find(abs(nodes(a,2))<1e-8);
    b3=find(abs(nodes(a,2)-(-c))<1e-8);
    bc([a(b1),a(b2),a(b3)],1) = true;
    bc(a(b2),2)=true;
    
    %% Reducción Matriz
    isFixed = reshape(bc',[],1);
    isFree = ~isFixed;
      
    % Solver
    Dr = K(isFree,isFree)\R(isFree);
    
    % Reconstrucción
    U = zeros(nDofTot,1);
    U(isFree) = U(isFree) + Dr;
    
    % Reacciones
    Rv = K(isFixed,isFree)*U(isFree);
    reacciones = nan(nDofTot,1);
    reacciones(isFixed) = Rv;
    reacciones = (reshape(reacciones,nDofNod,[]))';
    
    %% Recuperación de tensiones en los nodos
    stress = zeros(nel,nNodEle,3);
    uNod = [ -1 -1
        1 -1
        1  1
        -1  1 ];
    E_ene(m)=0;
    a   = 1/sqrt(3);
    % Ubicaciones puntos de Gauss
    upg = [ -a  -a
        a  -a
        a   a
        -a   a ];
    % Nï¿½mero de puntos de Gauss
    npg = size(upg,1);
    wpg = ones(npg,1);
    % Esquinas del elemento
    uNod = [ -1 -1
          1 -1
          1  1
         -1  1 ];
     
    E_sigma(m)=0;
    for iele = 1:nel
        nodesEle = nodes(elements(iele,:),:);
        for inode = 1:nNodEle
            % Punto de Gauss
            ksi = upg(ipg,1);
            eta = upg(ipg,2);
            % Derivadas de las funciones de forma respecto de ksi, eta
            dN = 1/4*[-(1-eta)   1-eta    1+eta  -(1+eta)
                -(1-ksi) -(1+ksi)   1+ksi    1-ksi ];
            % Derivadas de x,y, respecto de ksi, eta
            jac = dN*nodesEle;
            % Derivadas de las funciones de forma respecto de x,y.
            dNxy = jac\dN;          % dNxy = inv(jac)*dN
            
            B = zeros(size(C,2),nDofNod*nNodEle);
            B(1,1:2:7) = dNxy(1,:);
            B(2,2:2:8) = dNxy(2,:);
            B(3,1:2:7) = dNxy(2,:);
            B(3,2:2:8) = dNxy(1,:);
            
            eleDofs = nodeDofs(elements(iele,:),:);
            eleDofs = reshape(eleDofs',[],1);
            X0 = [mean(nodesEle(:,1)),mean(nodesEle(:,2))];
            % Ubico el punto de gauss local en la geometría espacial
            X = jac*[ksi,eta]' + X0';
            E_sigma(m) = E_sigma(m)+wpg(ipg)*(C*B*U(eleDofs) - ex_sol(X(1),X(2),'stress'))'*C_1*(C*B*U(eleDofs) - ex_sol(X(1),X(2),'stress'))*det(jac);
            E_ene(m) = E_ene(m)+wpg(ipg)*(B*U(eleDofs))'*C*B*U(eleDofs)*det(jac);
            
            % Para el Bandplot necesito los valores en las esquinas
            ksi = uNod(inode,1);
            eta = uNod(inode,2);
            % Derivadas de las funciones de forma respecto de ksi, eta
            dN = 1/4*[-(1-eta)   1-eta    1+eta  -(1+eta)
                -(1-ksi) -(1+ksi)   1+ksi    1-ksi ];
            % Derivadas de x,y, respecto de ksi, eta
            jac = dN*nodesEle;
            % Derivadas de las funciones de forma respecto de x,y.
            dNxy = jac\dN;          % dNxy = inv(jac)*dN
            B = zeros(size(C,2),nDofNod*nNodEle);
            B(1,1:2:7) = dNxy(1,:);
            B(2,2:2:8) = dNxy(2,:);
            B(3,1:2:7) = dNxy(2,:);
            B(3,2:2:8) = dNxy(1,:);
            stress(iele,inode,:) = C*B*U(eleDofs);
        end
    end
    %% Configuracion deformada
    U = (reshape(U,nDofNod,[]))';
    nodePosition = nodes + escala*U(:,1:2);
    Umax(m)=max(max(abs(U)));
    e_desp(m) = max(max(abs(Uex - U)))/max(max(abs(Uex)));
    % Gráficos
    % Sigmax
    figure()
    bandplot(elements,nodePosition,stress(:,:,1),[min(min(stress(:,:,1))),max(max(stress(:,:,1)))],'k');
    title('Sigma x')
    % Sigmay
%     figure()
%     bandplot(elements,nodePosition,stress(:,:,2),[min(min(stress(:,:,2))),max(max(stress(:,:,2)))],'k');
%     title('Sigma y')
    % txy
%     figure()
%     bandplot(elements,nodePosition,stress(:,:,3),[min(min(stress(:,:,3))),max(max(stress(:,:,3)))],'k');
%     title('Tau xy')
end
% Energía exacta
Eneex = integral2(@(x,y) ex_sol(x,y,'E'),0,L,-c,c);

txy = @(y) -(3*P)/(4*c)* (1-(y/c).^2);
sx = @(y) -(3*P*L*y)/(2*c^3);
% Fuerza total sobre borde izquierdo
Fyi = -integral(txy,-c,c);
Fxi = 0;

% Fuerza total sobre borde derecho
Fyd = integral(txy,-c,c);
Fxd = -integral(sx,-c,0)+-integral(sx,0,c);

figure
plot(h,E_ene,'-*r',h,Eneex*ones(length(h),1),'k')
xlabel('Parámetro de malla h')
ylabel('Energia')

figure
loglog(h,e_desp,'-*k')
xlabel('Parámetro de malla h')
ylabel('Error infinito relativo en desplazamientos')
[a,~] = polyfit(log(h(end-2:end)),log(e_desp(end-2:end)),1);
text((h(end)),(e_desp(end)),['Pendiente: ',num2str(a(1))]);

figure
loglog(h,sqrt(E_sigma),'-*k')
xlabel('Parámetro de malla h')
ylabel('Error en las tensiones')
[a,b] = polyfit(log(h(end-2:end)),log(sqrt(E_sigma(end-2:end))),1);
text((l(end)),(sqrt(E_sigma(end))),['Pendiente: ',num2str(a(1))]);