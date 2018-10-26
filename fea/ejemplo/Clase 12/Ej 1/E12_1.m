%% Práctica 12, ej 1
clear; close all; clc
% function feasolve(eleT,nodes,elements)
%% Datos
c = 10;
L = 100;
w = 1;
P = 80;
E = 1000;
nu = .25;
%eleT = 'Q4';

%% Nodos y elementos

nodes = nodes(:,2:3);
nDofNod = 2;                    % Número de grados de libertad por nodo
nNodEle = size(elements,2);     % Número de nodos por elemento
nel = size(elements,1);         % Número de elementos
nNod = size(nodes,1);           % Número de nodos
nDofTot = nDofNod*nNod;         % Número de grados de libertad

bc = false(nNod,nDofNod);       % Matriz de condiciones de borde
R = zeros(nNod,nDofNod);        % Vector de cargas

%meshplot(elements,nodes,'b')

%% Constitutiva
C = E/(1 - nu^2)*[ 1.0     nu         0.0
                   nu      1.0        0.0
                   0.0     0.0     (1 - nu)/2 ];
%% Gauss
switch eleT
    case 'Q4'
        [wpg, upg, npg] = gauss([2 2]);
    case 'Q8'
        [wpg, upg, npg] = gauss([3 3]);
    case 'Q9'
    case 'CST'
        a   = 1/2;
        upg = [a  0
               0  a
               a  a];    
        npg = size(upg,1);
        wpg = ones(npg,1)/3;
    case 'LST'
        a   = 1/2;
        upg = [a  0
               0  a
               a  a];    
        npg = size(upg,1);
        wpg = ones(npg,1)/3;
    case 'Q16'
end

%% Matriz de rigidez
K = zeros(nDofTot);
nodeDofs = reshape(1:nDofTot,nDofNod,nNod)';
for iele = 1:nel
    Ke = zeros(nDofNod*nNodEle);
    nodesEle = nodes(elements(iele,:),:);
    for ipg = 1:npg
        ksi = upg(ipg,1);
        eta = upg(ipg,2);  
        switch eleT
            case 'Q4'
                dN = shapefunsder([ksi eta],eleT);
            case 'Q8'
                dN = shapefunsder([ksi eta],eleT);
            case 'Q9'
            case 'CST'
                dN = [ -1, 1, 0; -1, 0, 1];
            case 'LST'
                dNLST =@(x,y) [ 4*x + 4*y - 3, 4*x - 1,       0, 4 - 4*y - 8*x, 4*y,          -4*y
                                4*x + 4*y - 3,       0, 4*y - 1,          -4*x, 4*x, 4 - 8*y - 4*x];
                dN = dNLST(ksi,eta);
            case 'Q16'
        end 
        jac = dN*nodesEle;                      
        dNxy = jac\dN;          % dNxy = inv(jac)*dN
        
        B = zeros(size(C,2),nDofNod*nNodEle);
        B(1,1:2:nNodEle*2-1) = dNxy(1,:);
        B(2,2:2:nNodEle*2) = dNxy(2,:);
        B(3,1:2:nNodEle*2-1) = dNxy(2,:);
        B(3,2:2:nNodEle*2) = dNxy(1,:); 

        Ke = Ke + B'*C*B*wpg(ipg)*det(jac);
    end
    eleDofs = nodeDofs(elements(iele,:),:);
    eleDofs = reshape(eleDofs',[],1);
    K(eleDofs,eleDofs) = K(eleDofs,eleDofs) + Ke;  
end

%% BC
bc(nodes(:,1)==100&nodes(:,2)==-10,1) =  true;
bc(nodes(:,1)==100&nodes(:,2)==0,:) =  true;
bc(nodes(:,1)==100&nodes(:,2)==10,1) =  true;

%% Cargas
sigx =@(x,y) -3/2*P*x*y/c^3;
tauxy =@(y) -3/4*P/c*(1-(y/c)^2);
npg = 2;
[wpg,upg] = gauss1D(npg);

switch eleT
    case 'Q4'
        for iele = 1:nel
            nodesEle = nodes(elements(iele,:),:);
            if nodesEle(2,1)==100
                for ipg = 1:npg
                    ksi = 1;
                    eta = upg(ipg);  
                    dN = shapefunsder([ksi eta],eleT);
                    N = shapefuns([ksi eta],eleT);
                    N = N([2 3]);
                    jac = dN*nodesEle;
                    sigma = [sigx(nodesEle(2,1),nodesEle(2,2)) 
                             sigx(nodesEle(3,1),nodesEle(3,2))];
                    tau = [tauxy(nodesEle(2,2)) 
                           tauxy(nodesEle(3,2))];
                    R(elements(iele,[2 3]),1) = R(elements(iele,[2 3]),1) + N'*N*w*wpg(ipg)*sigma*jac(2,2);
                    R(elements(iele,[2 3]),2) = R(elements(iele,[2 3]),2) + N'*N*w*wpg(ipg)*tau*jac(2,2);
                end
            end
            if nodesEle(1,1)==0
                for ipg = 1:npg
                    ksi = -1;
                    eta = upg(ipg);  
                    dN = shapefunsder([ksi eta],eleT);
                    N = shapefuns([ksi eta],eleT);
                    N = N([4 1]);
                    jac = dN*nodesEle;
                    tau = -[tauxy(nodesEle(4,2)) 
                            tauxy(nodesEle(1,2))];
                    R(elements(iele,[4 1]),2) = R(elements(iele,[4 1]),2) + N'*N*w*wpg(ipg)*tau*jac(2,2);
                end
            end
        end
    case 'Q8'
        for iele = 1:nel
            nodesEle = nodes(elements(iele,:),:);
            if nodesEle(2,1)==100
                for ipg = 1:npg
                    ksi = 1;
                    eta = upg(ipg);  
                    dN = shapefunsder([ksi eta],eleT);
                    N = shapefuns([ksi eta],eleT);
                    N = N([2 6 3]);
                    sigma = -[sigx(nodesEle(2,1),nodesEle(2,2)) 
                              sigx(nodesEle(6,1),nodesEle(6,2))
                              sigx(nodesEle(3,1),nodesEle(3,2))];
                    tau = [tauxy(nodesEle(2,2))
                           tauxy(nodesEle(6,2))
                           tauxy(nodesEle(3,2))];
                    jac = dN*nodesEle;
                    R(elements(iele,[2 6 3]),1) = R(elements(iele,[2 6 3]),1) + N'*N*w*wpg(ipg)*sigma*jac(2,2);
                    R(elements(iele,[2 6 3]),2) = R(elements(iele,[2 6 3]),2) + N'*N*w*wpg(ipg)*tau*jac(2,2);
                end
            end
            if nodesEle(1,1)==0
                for ipg = 1:npg
                    ksi = -1;
                    eta = upg(ipg);  
                    dN = shapefunsder([ksi eta],eleT);
                    N = shapefuns([ksi eta],eleT);
                    N = N([1 8 4]);
                    jac = dN*nodesEle;
                    tau = -[tauxy(nodesEle(1,2))
                            tauxy(nodesEle(8,2))
                            tauxy(nodesEle(4,2))];
                    R(elements(iele,[1 8 4]),2) = R(elements(iele,[1 8 4]),2) + N'*N*w*wpg(ipg)*tau*jac(2,2);
                end
            end
        end
    case 'Q9'
    case 'CST'
        for iele = 1:nel
            nodesEle = nodes(elements(iele,:),:);
            if nodesEle(3,1)==100
                for ipg = 1:npg
                    ksi = 1;
                    eta = upg(ipg);  
                    NCST =@(x,y) [ 1 - y - x, x, y];
                    N = NCST(ksi,eta);
                    dN = [ -1, 1, 0; -1, 0, 1];
                    N = N([2 3]);
                    sigma = -[sigx(nodesEle(2,1),nodesEle(2,2)) 
                              sigx(nodesEle(3,1),nodesEle(3,2))];
                    tau = [tauxy(nodesEle(2,2))
                           tauxy(nodesEle(3,2))];
                    jac = dN*nodesEle;
                    R(elements(iele,[2 3]),1) = R(elements(iele,[2 3]),1) + N'*N*w*wpg(ipg)*sigma*jac(2,2);
                    R(elements(iele,[2 3]),2) = R(elements(iele,[2 3]),2) + N'*N*w*wpg(ipg)*tau*jac(2,2);
                end
            end
             if nodesEle(3,1)==0
                for ipg = 1:npg
                    ksi = -1;
                    eta = upg(ipg);  
                    NCST =@(x,y) [ 1 - y - x, x, y];
                    N = NCST(ksi,eta);
                    dN = [ -1, 1, 0; -1, 0, 1];
                    N = N([1 3]);
                    tau = -[tauxy(nodesEle(1,2))
                            tauxy(nodesEle(3,2))];
                    jac = dN*nodesEle;
                    R(elements(iele,[1 3]),2) = R(elements(iele,[1 3]),2) + N'*N*w*wpg(ipg)*tau*jac(2,2);
                end
            end
        end
    case 'LST'
        for iele = 1:nel
            nodesEle = nodes(elements(iele,:),:);
            if nodesEle(3,1)==100
                for ipg = 1:npg
                    ksi = 1;
                    eta = upg(ipg);  
                    NLST =@(x,y) [ 2*x^2 + 4*x*y - 3*x + 2*y^2 - 3*y + 1, 2*x^2 - x, 2*y^2 - y, 4*x - 4*x*y - 4*x^2, 4*x*y, 4*y - 4*x*y - 4*y^2];
                    N = NLST(ksi,eta);
                    N = N([2 5 3]);
                    dNLST =@(x,y) [ 4*x + 4*y - 3, 4*x - 1,       0, 4 - 4*y - 8*x, 4*y,          -4*y
                                    4*x + 4*y - 3,       0, 4*y - 1,          -4*x, 4*x, 4 - 8*y - 4*x];
                    dN = dNLST(ksi,eta);
                    sigma = -[sigx(nodesEle(2,1),nodesEle(2,2)) 
                              sigx(nodesEle(5,1),nodesEle(5,2))
                              sigx(nodesEle(3,1),nodesEle(3,2))];
                    tau = [tauxy(nodesEle(2,2))
                           tauxy(nodesEle(5,2))
                           tauxy(nodesEle(3,2))];
                    jac = dN*nodesEle;
                    R(elements(iele,[2 5 3]),1) = R(elements(iele,[2 5 3]),1) + N'*N*w*wpg(ipg)*sigma*jac(2,2);
                    R(elements(iele,[2 5 3]),2) = R(elements(iele,[2 5 3]),2) + N'*N*w*wpg(ipg)*tau*jac(2,2);
                end
            end
            if nodesEle(3,1)==0
                for ipg = 1:npg
                    ksi = -1;
                    eta = upg(ipg);  
                    NLST =@(x,y) [ 2*x^2 + 4*x*y - 3*x + 2*y^2 - 3*y + 1, 2*x^2 - x, 2*y^2 - y, 4*x - 4*x*y - 4*x^2, 4*x*y, 4*y - 4*x*y - 4*y^2];
                    N = NLST(ksi,eta);
                    N = N([1 6 3]);
                    dNLST =@(x,y) [ 4*x + 4*y - 3, 4*x - 1,       0, 4 - 4*y - 8*x, 4*y,          -4*y
                                    4*x + 4*y - 3,       0, 4*y - 1,          -4*x, 4*x, 4 - 8*y - 4*x];
                    dN = dNLST(ksi,eta);
                    tau = -[tauxy(nodesEle(1,2))
                            tauxy(nodesEle(6,2))
                            tauxy(nodesEle(3,2))];
                    jac = dN*nodesEle;
                    R(elements(iele,[1 6 3]),2) = R(elements(iele,[1 6 3]),2) + N'*N*w*wpg(ipg)*tau*jac(2,2);
                end
            end
        end
    case 'Q16'
end

%% Reduccion Matriz
isFixed = reshape(bc',[],1);
isFree = ~isFixed;

Rr = reshape(R',[],1);

%% Solver
Dr = K(isFree,isFree)\Rr(isFree);

% Reconstrución
D = zeros(nDofTot,1);
D(isFree) = D(isFree) + Dr;

meshplot(elements,nodes+reshape(D,nDofNod,nNod)','r')
Dr = reshape(D,nDofNod,nNod)';

%% Deformación en los nodos
epsilon = zeros(nel,nNodEle,3);
switch eleT
    case 'Q4'
        unod = [-1 -1
                1 -1
                1 1
                -1 1];
end
for iele = 1:nel
    nodesEle = nodes(elements(iele,:),:);
    for in = 1:nNodEle
        ksi = unod(in,1);
        eta = unod(in,2);  
        switch eleT
            case 'Q4'
                dN = shapefunsder([ksi eta],eleT);
            case 'Q8'
                dN = shapefunsder([ksi eta],eleT);
            case 'Q9'
            case 'CST'
                dN = [ -1, 1, 0; -1, 0, 1];
            case 'LST'
                dNLST =@(x,y) [ 4*x + 4*y - 3, 4*x - 1,       0, 4 - 4*y - 8*x, 4*y,          -4*y
                                4*x + 4*y - 3,       0, 4*y - 1,          -4*x, 4*x, 4 - 8*y - 4*x];
                dN = dNLST(ksi,eta);
            case 'Q16'
        end 
        jac = dN*nodesEle;                      
        dNxy = jac\dN;
        
        B = zeros(size(C,2),nDofNod*nNodEle);
        B(1,1:2:nNodEle*2-1) = dNxy(1,:);
        B(2,2:2:nNodEle*2) = dNxy(2,:);
        B(3,1:2:nNodEle*2-1) = dNxy(2,:);
        B(3,2:2:nNodEle*2) = dNxy(1,:); 
        eleDofs = nodeDofs(elements(iele,:),:);
        eleDofs = reshape(eleDofs',[],1);
        epsilon(iele,in,:) = B*D(eleDofs);
    end
end

%% Solución exacta

syms x y real
I = (2*c)^3*w/12;
G = E/2/(1+nu);
u = P/(2*I)*(-(x^2-L^2)/E*y-nu*y*(y^2-c^2)/(3*E)+y*(y^2-c^2)/(3*G));
v = P/I*(nu*x*y^2/(2*E)+(x^3-L^3)/(6*E)-(L^2/(2*E)+nu*c^2/(6*E)+c^2/(3*G))*(x-L));
e11 = diff(u,x);
e22 = diff(v,y);
e12 = (diff(u,y)+diff(v,x))/2;

u = matlabFunction(u);
v = matlabFunction(v);
e11 = matlabFunction(e11);
e22 = matlabFunction(e22);
e12 = matlabFunction(e12);

Dexact = [u(nodes(:,1),nodes(:,2)) v(nodes(:,1),nodes(:,2))];
hold on
meshplot(elements,nodes+Dexact,'g')

epsilonExact = zeros(nel,nNodEle,3);
for iele = 1:nel
    nodesEle = nodes(elements(iele,:),:);
    for in = 1:nNodEle
        epsilonExact(iele,in,:) = [e11(nodesEle(in,1),nodesEle(in,2)) ...
            e22(nodesEle(in,1),nodesEle(in,2)) e12(nodesEle(in,2))];
    end
end

%% Gauss
switch eleT
    case 'Q4'
        [wpg, upg, npg] = gauss([2 2]);
    case 'Q8'
        [wpg, upg, npg] = gauss([3 3]);
    case 'Q9'
    case 'CST'
        a   = 1/2;
        upg = [a  0
               0  a
               a  a];    
        npg = size(upg,1);
        wpg = ones(npg,1)/3;
    case 'LST'
        a   = 1/2;
        upg = [a  0
               0  a
               a  a];    
        npg = size(upg,1);
        wpg = ones(npg,1)/3;
    case 'Q16'
end

%% Error
Error = 0;
delta = epsilon-epsilonExact;
for iele = 1:nel
    nodesEle = nodes(elements(iele,:),:);
    d = zeros(1,3);
    for ipg = 1:npg
        ksi = upg(ipg,1);
        eta = upg(ipg,2);  
        switch eleT
            case 'Q4'
                dN = shapefunsder([ksi eta],eleT);
            case 'Q8'
                dN = shapefunsder([ksi eta],eleT);
            case 'Q9'
            case 'CST'
                dN = [ -1, 1, 0; -1, 0, 1];
            case 'LST'
                dNLST =@(x,y) [ 4*x + 4*y - 3, 4*x - 1,       0, 4 - 4*y - 8*x, 4*y,          -4*y
                                4*x + 4*y - 3,       0, 4*y - 1,          -4*x, 4*x, 4 - 8*y - 4*x];
                dN = dNLST(ksi,eta);
            case 'Q16'
        end 
        jac = dN*nodesEle;                      
        dNxy = jac\dN;          % dNxy = inv(jac)*dN
        dpg = [delta(iele,ipg,1) delta(iele,ipg,2) delta(iele,ipg,3)];
        d = d + dpg*w*wpg(ipg)*det(jac);
    end
    Error = Error + d*C*d';
end












% U = 0;
% for iele = 1:nel
%     nodesEle = nodes(elements(iele,:),:);
%     d = zeros(1,3);
%     for ipg = 1:npg
%         ksi = upg(ipg,1);
%         eta = upg(ipg,2);  
%         switch eleT
%             case 'Q4'
%                 dN = shapefunsder([ksi eta],eleT);
%             case 'Q8'
%                 dN = shapefunsder([ksi eta],eleT);
%             case 'Q9'
%             case 'CST'
%                 dN = [ -1, 1, 0; -1, 0, 1];
%             case 'LST'
%                 dNLST =@(x,y) [ 4*x + 4*y - 3, 4*x - 1,       0, 4 - 4*y - 8*x, 4*y,          -4*y
%                                 4*x + 4*y - 3,       0, 4*y - 1,          -4*x, 4*x, 4 - 8*y - 4*x];
%                 dN = dNLST(ksi,eta);
%             case 'Q16'
%         end 
%         jac = dN*nodesEle;                      
%         dNxy = jac\dN;          % dNxy = inv(jac)*dN
%         dpg = [epsilon(iele,ipg,1) epsilon(iele,ipg,2) epsilon(iele,ipg,3)];
%         d = d + dpg*w*wpg(ipg)*det(jac);
%     end
%     U = U + d*C*d';
% end





% end






