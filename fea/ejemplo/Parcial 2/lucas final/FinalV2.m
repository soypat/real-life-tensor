%% Inicializacion
tic
clc
clear
Caso='Normal';
Direccion='x';

E1=200000;
NU1=0.27;
E2=117000;
NU2=0.34;
Espesor1=0.4;

        alfa1=0.000011;
        T=100;
        alfa2=1.65e-005;
        LadoH=7;
        Largo=30;
        Espesor=2;
        ElementosH=14;
        ElementosLargo=60;
        ElementosEspesor=10;


%% PatchTest y Mallador
switch Caso
    %% Patch Test
    
    case 'Patch'
        T=0;
        Largo=2;
        Espesor=2;
        LadoH=2;
        ElementosH=2;
        ElementosLargo=2;
        ElementosEspesor=2;
        Espesor1=2;
        format short g
        
        
        H=0:LadoH/ElementosH:LadoH;
        L=0:Largo/ElementosLargo:Largo;
        W=0:Espesor/ElementosEspesor:Espesor;
        NH=length(H);
        NL=length(L);
        NW=length(W);
        nodes=zeros(NL*NH*NW,3);
        for j=1:length(W)
            for i=1:length(H)
                nodes(NH*NL*(j-1)+NL*(i-1)+[1:NL],:)=nodes(NH*NL*(j-1)+NL*(i-1)+[1:NL],:)+[L' ones(NL,1)*H(i) ones(NL,1)*W(j)];
            end
        end
        elements=zeros(ElementosH*ElementosLargo*ElementosEspesor,8);
        erase=0;
        for j=1:length(W)-1
            elements(((ElementosH*(ElementosLargo+1)*(j-1)+1):(ElementosH*(ElementosLargo+1)*j-1)),:)=[(1:(NL*(NH-1)-1))+NH*NL*(j-1) ; (2:(NL*(NH-1)))+NH*NL*(j-1) ; ((NL+2):(NL*NH))+NH*NL*(j-1) ; ((NL+1):((NL*NH)-1))+NH*NL*(j-1);(1:(NL*(NH-1)-1))+NH*NL*(j); (2:(NL*(NH-1)))+NH*NL*(j) ; ((NL+2):(NL*NH))+NH*NL*(j) ; ((NL+1):((NL*NH)-1))+NH*NL*(j)]';
            erase=[erase [NL:NL:NL*(NH-1)]+(j-1)*ElementosH*(ElementosLargo+1)];
        end
        erase(end)=[];
        erase(erase==0)=[];
        elements(erase,:)=[];
        for i=1:ElementosH*ElementosLargo*ElementosEspesor
            node1=nodes(elements(i,1),:);
            node2=nodes(elements(i,2),:);
            node3=nodes(elements(i,3),:);
            node4=nodes(elements(i,4),:);
            node5=nodes(elements(i,5),:);
            node6=nodes(elements(i,6),:);
            node7=nodes(elements(i,7),:);
            node8=nodes(elements(i,8),:);
            nodes=[nodes; node1+(node2-node1)./2;node2+(node3-node2)/2;node3+(node4-node3)/2;node4+(node1-node4)/2;node1+(node5-node1)/2;node2+(node6-node2)/2;node3+(node7-node3)/2;node4+(node8-node4)/2;node5+(node6-node5)/2;node6+(node7-node6)/2;node7+(node8-node7)/2;node8+(node5-node8)/2];
        end
        nodes=sortrows(unique(nodes,'rows'),[3 2 1]);
        H20elements=zeros(ElementosLargo*ElementosH*ElementosEspesor,20);
        for j=1:ElementosEspesor
            for i=1:ElementosH
                Pos=(ElementosLargo*(i-1)+(ElementosLargo*ElementosH)*(j-1)+1):(ElementosLargo*(i)+(ElementosLargo*ElementosH)*(j-1));
                H20elements(Pos,1)=[((1:2:2*NL-2)+(i-1)*((NL*2-1)+NL))+(j-1)*(NH*NL+((2*NL-1)*(NH))+ElementosH*NL)]';
                H20elements(Pos,4)=[(1:2:2*NL-2)+(i)*((NL*2-1)+NL)+(j-1)*(NH*NL+((2*NL-1)*(NH))+ElementosH*NL)]';
                H20elements(Pos,5)=[((1:2:2*NL-2)+(i-1)*((NL*2-1)+NL))+(j)*(NH*NL+((2*NL-1)*(NH))+ElementosH*NL)]';
                H20elements(Pos,8)=[(1:2:2*NL-2)+(i)*((NL*2-1)+NL)+(j)*(NH*NL+((2*NL-1)*(NH))+ElementosH*NL)]';
                H20elements(Pos,10)=[(2*NL+1:3*NL-1)+(i-1)*((NL*2-1)+NL)+(j-1)*(NH*NL+((2*NL-1)*(NH))+ElementosH*NL)]';
                H20elements(Pos,13)=[((2*NL+(NH-1)*(NL*3-1)):2*NL+(NH-1)*(NL*3-1)+NL-2)+NL*(i-1)+(j-1)*(NH*NL-1+(2*NL+(NH-1)*(NL*3-1)))]';
                H20elements(Pos,18)=[(2*NL+1:3*NL-1)+(i-1)*((NL*2-1)+NL)+(j)*(NH*NL+((2*NL-1)*(NH))+ElementosH*NL)]';
                
            end
        end
        H20elements(:,2)=H20elements(:,1)+2;
        H20elements(:,3)=H20elements(:,4)+2;
        H20elements(:,6)=H20elements(:,5)+2;
        H20elements(:,7)=H20elements(:,8)+2;
        H20elements(:,9)=H20elements(:,1)+1;
        H20elements(:,11)=H20elements(:,4)+1;
        H20elements(:,12)=H20elements(:,10)-1;
        H20elements(:,14)=H20elements(:,13)+1;
        H20elements(:,16)=H20elements(:,13)+NL;
        H20elements(:,15)=H20elements(:,16)+1;
        H20elements(:,17)=H20elements(:,5)+1;
        H20elements(:,19)=H20elements(:,8)+1;
        H20elements(:,20)=H20elements(:,18)-1;
        elements=H20elements;
        nodes(floor(length(nodes)/2)+1,:)=[0.75 1.15 1.36];
        nDofNod = 3;                    % Nï¿½mero de grados de libertad por nodo
        nNodEle = 20;                   % Nï¿½mero de nodos por elemento
        nel = size(elements,1);         % Nï¿½mero de elementos
        nNod = size(nodes,1);           % Nï¿½mero de nodos
        nDofTot = nDofNod*nNod;         % Nï¿½mero de grados de libertad
        R = zeros(nNod,nDofNod);        % Matriz de cargas
        bc = zeros(nNod,nDofNod);        % Matriz de condiciones de borde
        
        bc(1,:)=1;
        bc(5,2:3)=1;
        bc(6,3)=1;
        switch Direccion
            case 'z'
                for iele = 1:4
                    
                    R(elements(iele,1:4),3)=R(elements(iele,1:4),3)-1/12;
                    
                    R(elements(iele,(1:4)+8),3)=R(elements(iele,(1:4)+8),3)+1/3;
                    
                    R(elements(iele+4,(1:4)+4),3)=R(elements(iele+4,(1:4)+4),3)+1/12;
                    
                    R(elements(iele+4,(1:4)+16),3)=R(elements(iele+4,(1:4)+16),3)-1/3;
                end
                
            case 'x'
                for iele = 1:2:7
                    
                    R(elements(iele,[1 4 5 8]),1)=R(elements(iele,[1 4 5 8]),1)-1/12;
                    
                    R(elements(iele,[12 13 16 20]),1)=R(elements(iele,[12 13 16 20]),1)+1/3;
                    
                    R(elements(iele+1,[2 3 6 7]),1)=R(elements(iele+1,[2 3 6 7]),1)+1/12;
                    
                    R(elements(iele+1,[10 14 15 18]),1)=R(elements(iele+1,[10 14 15 18]),1)-1/3;
                end
            case 'y'
                for iele = [1 2 5 6]
                    
                    R(elements(iele,[1 2 5 6]),2)=R(elements(iele,[1 2 5 6]),2)-1/12;
                    
                    R(elements(iele,[9 13 14 17]),2)=R(elements(iele,[9 13 14 17]),2)+1/3;
                    
                    R(elements(iele+2,[3 4 7 8]),2)=R(elements(iele+2,[3 4 7 8]),2)+1/12;
                    
                    R(elements(iele+2,[11 15 16 19]),2)=R(elements(iele+2,[11 15 16 19]),2)-1/3;
                end
            case 'yz'
                for iele = [1 2 5 6]
                    
                    R(elements(iele,[1 2 5 6]),3)=R(elements(iele,[1 2 5 6]),3)-1/12;
                    
                    R(elements(iele,[9 13 14 17]),3)=R(elements(iele,[9 13 14 17]),3)+1/3;
                    
                    R(elements(iele+2,[3 4 7 8]),3)=R(elements(iele+2,[3 4 7 8]),3)+1/12;
                    
                    R(elements(iele+2,[11 15 16 19]),3)=R(elements(iele+2,[11 15 16 19]),3)-1/3;
                end
                bc = zeros(nNod,nDofNod);
                bc(nodes(:,3)==0,:)=1;
                for iele = 1:4
                    
                    R(elements(iele,1:4),2)=R(elements(iele,1:4),2)-1/12;
                    
                    R(elements(iele,(1:4)+8),2)=R(elements(iele,(1:4)+8),2)+1/3;
                    
                    R(elements(iele+4,(1:4)+4),2)=R(elements(iele+4,(1:4)+4),2)+1/12;
                    
                    R(elements(iele+4,(1:4)+16),2)=R(elements(iele+4,(1:4)+16),2)-1/3;
                end
                R=-R;
            case 'xz'
                for iele = 1:2:7
                    
                    R(elements(iele,[1 4 5 8]),3)=R(elements(iele,[1 4 5 8]),3)-1/12;
                    
                    R(elements(iele,[12 13 16 20]),3)=R(elements(iele,[12 13 16 20]),3)+1/3;
                    
                    R(elements(iele+1,[2 3 6 7]),3)=R(elements(iele+1,[2 3 6 7]),3)+1/12;
                    
                    R(elements(iele+1,[10 14 15 18]),3)=R(elements(iele+1,[10 14 15 18]),3)-1/3;
                end
                bc = zeros(nNod,nDofNod);
                bc(nodes(:,3)==0,:)=1;
                for iele = 1:4
                    
                    R(elements(iele,1:4),1)=R(elements(iele,1:4),1)-1/12;
                    
                    R(elements(iele,(1:4)+8),1)=R(elements(iele,(1:4)+8),1)+1/3;
                    
                    R(elements(iele+4,(1:4)+4),1)=R(elements(iele+4,(1:4)+4),1)+1/12;
                    
                    R(elements(iele+4,(1:4)+16),1)=R(elements(iele+4,(1:4)+16),1)-1/3;
                end
                R=-R;
            case 'xy'
                for iele = 1:2:7
                    
                    R(elements(iele,[1 4 5 8]),2)=R(elements(iele,[1 4 5 8]),2)-1/12;
                    
                    R(elements(iele,[12 13 16 20]),2)=R(elements(iele,[12 13 16 20]),2)+1/3;
                    
                    R(elements(iele+1,[2 3 6 7]),2)=R(elements(iele+1,[2 3 6 7]),2)+1/12;
                    
                    R(elements(iele+1,[10 14 15 18]),2)=R(elements(iele+1,[10 14 15 18]),2)-1/3;
                end
                bc = zeros(nNod,nDofNod);
                bc(nodes(:,2)==0,:)=1;
                for iele = [1 2 5 6]
                    
                    R(elements(iele,[1 2 5 6]),1)=R(elements(iele,[1 2 5 6]),1)-1/12;
                    
                    R(elements(iele,[9 13 14 17]),1)=R(elements(iele,[9 13 14 17]),1)+1/3;
                    
                    R(elements(iele+2,[3 4 7 8]),1)=R(elements(iele+2,[3 4 7 8]),1)+1/12;
                    
                    R(elements(iele+2,[11 15 16 19]),1)=R(elements(iele+2,[11 15 16 19]),1)-1/3;
                end
                R=-R;
        end
        
        
%          R=R*12;
%         scatter3(nodes(:,1),nodes(:,2),nodes(:,3))
%         hold on
%         meshplot(elements,nodes,'r')
        %% Bi-metallic Bar
    case 'Normal'
        
        

 
        
        
        H=0:LadoH/ElementosH:LadoH;
        L=0:Largo/ElementosLargo:Largo;
        W=0:Espesor/ElementosEspesor:Espesor;
        NH=length(H);
        NL=length(L);
        NW=length(W);
        nodes=zeros(NL*NH*NW,3);
        for j=1:length(W)
            for i=1:length(H)
                nodes(NH*NL*(j-1)+NL*(i-1)+[1:NL],:)=nodes(NH*NL*(j-1)+NL*(i-1)+[1:NL],:)+[L' ones(NL,1)*H(i) ones(NL,1)*W(j)];
            end
        end
        elements=zeros(ElementosH*ElementosLargo*ElementosEspesor,8);
        erase=0;
        for j=1:length(W)-1
            elements(((ElementosH*(ElementosLargo+1)*(j-1)+1):(ElementosH*(ElementosLargo+1)*j-1)),:)=[(1:(NL*(NH-1)-1))+NH*NL*(j-1) ; (2:(NL*(NH-1)))+NH*NL*(j-1) ; ((NL+2):(NL*NH))+NH*NL*(j-1) ; ((NL+1):((NL*NH)-1))+NH*NL*(j-1);(1:(NL*(NH-1)-1))+NH*NL*(j); (2:(NL*(NH-1)))+NH*NL*(j) ; ((NL+2):(NL*NH))+NH*NL*(j) ; ((NL+1):((NL*NH)-1))+NH*NL*(j)]';
            erase=[erase [NL:NL:NL*(NH-1)]+(j-1)*ElementosH*(ElementosLargo+1)];
        end
        erase(end)=[];
        erase(erase==0)=[];
        elements(erase,:)=[];
        for i=1:ElementosH*ElementosLargo*ElementosEspesor
            node1=nodes(elements(i,1),:);
            node2=nodes(elements(i,2),:);
            node3=nodes(elements(i,3),:);
            node4=nodes(elements(i,4),:);
            node5=nodes(elements(i,5),:);
            node6=nodes(elements(i,6),:);
            node7=nodes(elements(i,7),:);
            node8=nodes(elements(i,8),:);
            nodes=[nodes; node1+(node2-node1)./2;node2+(node3-node2)/2;node3+(node4-node3)/2;node4+(node1-node4)/2;node1+(node5-node1)/2;node2+(node6-node2)/2;node3+(node7-node3)/2;node4+(node8-node4)/2;node5+(node6-node5)/2;node6+(node7-node6)/2;node7+(node8-node7)/2;node8+(node5-node8)/2];
        end
        Nodes=sortrows(unique(nodes,'rows'),[3 2 1]);
        H20elements=zeros(ElementosLargo*ElementosH*ElementosEspesor,20);
        for j=1:ElementosEspesor
            for i=1:ElementosH
                Pos=(ElementosLargo*(i-1)+(ElementosLargo*ElementosH)*(j-1)+1):(ElementosLargo*(i)+(ElementosLargo*ElementosH)*(j-1));
                H20elements(Pos,1)=[((1:2:2*NL-2)+(i-1)*((NL*2-1)+NL))+(j-1)*(NH*NL+((2*NL-1)*(NH))+ElementosH*NL)]';
                H20elements(Pos,4)=[(1:2:2*NL-2)+(i)*((NL*2-1)+NL)+(j-1)*(NH*NL+((2*NL-1)*(NH))+ElementosH*NL)]';
                H20elements(Pos,5)=[((1:2:2*NL-2)+(i-1)*((NL*2-1)+NL))+(j)*(NH*NL+((2*NL-1)*(NH))+ElementosH*NL)]';
                H20elements(Pos,8)=[(1:2:2*NL-2)+(i)*((NL*2-1)+NL)+(j)*(NH*NL+((2*NL-1)*(NH))+ElementosH*NL)]';
                H20elements(Pos,10)=[(2*NL+1:3*NL-1)+(i-1)*((NL*2-1)+NL)+(j-1)*(NH*NL+((2*NL-1)*(NH))+ElementosH*NL)]';
                H20elements(Pos,13)=[((2*NL+(NH-1)*(NL*3-1)):2*NL+(NH-1)*(NL*3-1)+NL-2)+NL*(i-1)+(j-1)*(NH*NL-1+(2*NL+(NH-1)*(NL*3-1)))]';
                H20elements(Pos,18)=[(2*NL+1:3*NL-1)+(i-1)*((NL*2-1)+NL)+(j)*(NH*NL+((2*NL-1)*(NH))+ElementosH*NL)]';
                
            end
        end
        H20elements(:,2)=H20elements(:,1)+2;
        H20elements(:,3)=H20elements(:,4)+2;
        H20elements(:,6)=H20elements(:,5)+2;
        H20elements(:,7)=H20elements(:,8)+2;
        H20elements(:,9)=H20elements(:,1)+1;
        H20elements(:,11)=H20elements(:,4)+1;
        H20elements(:,12)=H20elements(:,10)-1;
        H20elements(:,14)=H20elements(:,13)+1;
        H20elements(:,16)=H20elements(:,13)+NL;
        H20elements(:,15)=H20elements(:,16)+1;
        H20elements(:,17)=H20elements(:,5)+1;
        H20elements(:,19)=H20elements(:,8)+1;
        H20elements(:,20)=H20elements(:,18)-1;
        elements=H20elements;
        nodes=Nodes;
        %%
        nDofNod = 3;                    % Nï¿½mero de grados de libertad por nodo
        nNodEle = 20;                   % Nï¿½mero de nodos por elemento
        nel = size(elements,1);         % Nï¿½mero de elementos
        nNod = size(nodes,1);           % Nï¿½mero de nodos
        nDofTot = nDofNod*nNod;         % Nï¿½mero de grados de libertad
        R = zeros(nNod,nDofNod);        % Matriz de cargas
        bc=zeros(nNod,3);
        bc(find(nodes(:,1)==0),:)=1;
        % bc(5,2:3)=1;
        % bc(6,3)=1;
        %
        
end





%% Matriz Constitutiva Metal 1
Lambda1=(E1*NU1)/((1+NU1)*(1-2*NU1));
MU1=E1/(2+2*NU1);
C1 = [Lambda1+2*MU1   Lambda1         Lambda1       0   0   0
    Lambda1         Lambda1+2*MU1   Lambda1       0   0   0
    Lambda1         Lambda1         Lambda1+2*MU1 0   0   0
    0               0               0             MU1 0   0
    0               0               0             0   MU1 0
    0               0               0             0   0   MU1];

% C1=sparse(C1);



%% Matriz Constitutiva Metal 2
Lambda2=(E2*NU2)/((1+NU2)*(1-2*NU2));
MU2=E2/(2+2*NU2);
C2 = [Lambda2+2*MU2     Lambda2          Lambda2        0   0   0
    Lambda2           Lambda2+2*MU2    Lambda2        0   0   0
    Lambda2           Lambda2          Lambda2+2*MU2  0   0   0
    0                 0                0              MU2 0   0
    0                 0                0              0   MU2 0
    0                 0                0              0   0   MU2];

% C2=sparse(C2);



%% Puntos de Gauss Full
a=(3/5)^0.5;
c1=5/9;
c2=8/9;
xyzgauss=[-1    -1   -1
        1    -1   -1
        1     1   -1
       -1     1   -1
       -1    -1    1
        1    -1    1
        1     1    1
       -1     1    1
        0    -1   -1
        1     0   -1
        0     1   -1
       -1     0   -1
       -1    -1    0
        1    -1    0
        1     1    0
       -1     1    0
        0    -1    1
        1     0    1
        0     1    1
       -1     0    1
        0     0   -1
       -1     0    0
        0    -1    0
        1     0    0
        0     1    0
        0     0    1
        0     0    0];
upg=a*[xyzgauss]; % Ubicaciones puntos de Gauss
wpg=zeros(27,1);
wpg([1 2 3 4 5 6 7 8])=c1^3;
wpg([9 10 11 12 13 14 15 16 17 18 19 20])=c1^2*c2;
wpg([21 22 23 24 25 26])=c1*c2^2;
wpg(27)=c2^3; % Pesos puntos de Gauss
npg = size(upg,1); % Numero de puntos de Gauss




%% Matriz de rigidez
K = spalloc(nDofTot,nDofTot,nDofTot*180);
F = spalloc(nDofTot,1,ceil(nDofTot/2));
nodeDofs = reshape(1:nDofTot,nDofNod,nNod)';
count=0;
for iele = 1:nel
    Fe=spalloc(60,1,60);
    Ke = zeros(nDofNod*nNodEle,nDofNod*nNodEle);
    nodesEle = nodes(elements(iele,:),:);
    if max(nodes(elements(iele,:),3))<= Espesor1
        C=C1;
        epsilon=alfa1*T*[1;1;1;0;0;0];
    else
        C=C2;
        epsilon=alfa2*T*[1;1;1;0;0;0];
    end
    for ipg = 1:npg
        count=count+1;
        % Punto de Gauss
        ksi  = upg(ipg,1);
        eta  = upg(ipg,2);
        zeta = upg(ipg,3);
        % Derivadas de las funciones de forma respecto de ksi, eta
        % [N dN]=FuncionFormaH20(ksi,eta,zeta);
        N =[ (eta^2*ksi*zeta)/8 - (eta^2*ksi)/8 - (eta^2*zeta)/8 + eta^2/8 + (eta*ksi^2*zeta)/8 - (eta*ksi^2)/8 + (eta*ksi*zeta^2)/8 - (eta*ksi*zeta)/8 - (eta*zeta^2)/8 + eta/8 - (ksi^2*zeta)/8 + ksi^2/8 - (ksi*zeta^2)/8 + ksi/8 + zeta^2/8 + zeta/8 - 1/4, - (eta^2*ksi*zeta)/8 + (eta^2*ksi)/8 - (eta^2*zeta)/8 + eta^2/8 + (eta*ksi^2*zeta)/8 - (eta*ksi^2)/8 - (eta*ksi*zeta^2)/8 + (eta*ksi*zeta)/8 - (eta*zeta^2)/8 + eta/8 - (ksi^2*zeta)/8 + ksi^2/8 + (ksi*zeta^2)/8 - ksi/8 + zeta^2/8 + zeta/8 - 1/4, - (eta^2*ksi*zeta)/8 + (eta^2*ksi)/8 - (eta^2*zeta)/8 + eta^2/8 - (eta*ksi^2*zeta)/8 + (eta*ksi^2)/8 + (eta*ksi*zeta^2)/8 - (eta*ksi*zeta)/8 + (eta*zeta^2)/8 - eta/8 - (ksi^2*zeta)/8 + ksi^2/8 + (ksi*zeta^2)/8 - ksi/8 + zeta^2/8 + zeta/8 - 1/4, (eta^2*ksi*zeta)/8 - (eta^2*ksi)/8 - (eta^2*zeta)/8 + eta^2/8 - (eta*ksi^2*zeta)/8 + (eta*ksi^2)/8 - (eta*ksi*zeta^2)/8 + (eta*ksi*zeta)/8 + (eta*zeta^2)/8 - eta/8 - (ksi^2*zeta)/8 + ksi^2/8 - (ksi*zeta^2)/8 + ksi/8 + zeta^2/8 + zeta/8 - 1/4, - (eta^2*ksi*zeta)/8 - (eta^2*ksi)/8 + (eta^2*zeta)/8 + eta^2/8 - (eta*ksi^2*zeta)/8 - (eta*ksi^2)/8 + (eta*ksi*zeta^2)/8 + (eta*ksi*zeta)/8 - (eta*zeta^2)/8 + eta/8 + (ksi^2*zeta)/8 + ksi^2/8 - (ksi*zeta^2)/8 + ksi/8 + zeta^2/8 - zeta/8 - 1/4, (eta^2*ksi*zeta)/8 + (eta^2*ksi)/8 + (eta^2*zeta)/8 + eta^2/8 - (eta*ksi^2*zeta)/8 - (eta*ksi^2)/8 - (eta*ksi*zeta^2)/8 - (eta*ksi*zeta)/8 - (eta*zeta^2)/8 + eta/8 + (ksi^2*zeta)/8 + ksi^2/8 + (ksi*zeta^2)/8 - ksi/8 + zeta^2/8 - zeta/8 - 1/4, (eta^2*ksi*zeta)/8 + (eta^2*ksi)/8 + (eta^2*zeta)/8 + eta^2/8 + (eta*ksi^2*zeta)/8 + (eta*ksi^2)/8 + (eta*ksi*zeta^2)/8 + (eta*ksi*zeta)/8 + (eta*zeta^2)/8 - eta/8 + (ksi^2*zeta)/8 + ksi^2/8 + (ksi*zeta^2)/8 - ksi/8 + zeta^2/8 - zeta/8 - 1/4, - (eta^2*ksi*zeta)/8 - (eta^2*ksi)/8 + (eta^2*zeta)/8 + eta^2/8 + (eta*ksi^2*zeta)/8 + (eta*ksi^2)/8 - (eta*ksi*zeta^2)/8 - (eta*ksi*zeta)/8 + (eta*zeta^2)/8 - eta/8 + (ksi^2*zeta)/8 + ksi^2/8 - (ksi*zeta^2)/8 + ksi/8 + zeta^2/8 - zeta/8 - 1/4, (eta*zeta)/4 - zeta/4 - eta/4 + (eta*ksi^2)/4 + (ksi^2*zeta)/4 - ksi^2/4 - (eta*ksi^2*zeta)/4 + 1/4, ksi/4 - zeta/4 - (ksi*zeta)/4 - (eta^2*ksi)/4 + (eta^2*zeta)/4 - eta^2/4 + (eta^2*ksi*zeta)/4 + 1/4, eta/4 - zeta/4 - (eta*zeta)/4 - (eta*ksi^2)/4 + (ksi^2*zeta)/4 - ksi^2/4 + (eta*ksi^2*zeta)/4 + 1/4, (ksi*zeta)/4 - zeta/4 - ksi/4 + (eta^2*ksi)/4 + (eta^2*zeta)/4 - eta^2/4 - (eta^2*ksi*zeta)/4 + 1/4, (eta*ksi)/4 - ksi/4 - eta/4 + (eta*zeta^2)/4 + (ksi*zeta^2)/4 - zeta^2/4 - (eta*ksi*zeta^2)/4 + 1/4, ksi/4 - eta/4 - (eta*ksi)/4 + (eta*zeta^2)/4 - (ksi*zeta^2)/4 - zeta^2/4 + (eta*ksi*zeta^2)/4 + 1/4, eta/4 + ksi/4 + (eta*ksi)/4 - (eta*zeta^2)/4 - (ksi*zeta^2)/4 - zeta^2/4 - (eta*ksi*zeta^2)/4 + 1/4, eta/4 - ksi/4 - (eta*ksi)/4 - (eta*zeta^2)/4 + (ksi*zeta^2)/4 - zeta^2/4 + (eta*ksi*zeta^2)/4 + 1/4, zeta/4 - eta/4 - (eta*zeta)/4 + (eta*ksi^2)/4 - (ksi^2*zeta)/4 - ksi^2/4 + (eta*ksi^2*zeta)/4 + 1/4, ksi/4 + zeta/4 + (ksi*zeta)/4 - (eta^2*ksi)/4 - (eta^2*zeta)/4 - eta^2/4 - (eta^2*ksi*zeta)/4 + 1/4, eta/4 + zeta/4 + (eta*zeta)/4 - (eta*ksi^2)/4 - (ksi^2*zeta)/4 - ksi^2/4 - (eta*ksi^2*zeta)/4 + 1/4, zeta/4 - ksi/4 - (ksi*zeta)/4 + (eta^2*ksi)/4 - (eta^2*zeta)/4 - eta^2/4 + (eta^2*ksi*zeta)/4 + 1/4];
        dN=[ ksi/4 - (eta*ksi)/4 - (eta*zeta)/8 - (ksi*zeta)/4 + (eta*zeta^2)/8 + (eta^2*zeta)/8 - eta^2/8 - zeta^2/8 + (eta*ksi*zeta)/4 + 1/8, ksi/4 - (eta*ksi)/4 + (eta*zeta)/8 - (ksi*zeta)/4 - (eta*zeta^2)/8 - (eta^2*zeta)/8 + eta^2/8 + zeta^2/8 + (eta*ksi*zeta)/4 - 1/8, ksi/4 + (eta*ksi)/4 - (eta*zeta)/8 - (ksi*zeta)/4 + (eta*zeta^2)/8 - (eta^2*zeta)/8 + eta^2/8 + zeta^2/8 - (eta*ksi*zeta)/4 - 1/8, ksi/4 + (eta*ksi)/4 + (eta*zeta)/8 - (ksi*zeta)/4 - (eta*zeta^2)/8 + (eta^2*zeta)/8 - eta^2/8 - zeta^2/8 - (eta*ksi*zeta)/4 + 1/8, ksi/4 - (eta*ksi)/4 + (eta*zeta)/8 + (ksi*zeta)/4 + (eta*zeta^2)/8 - (eta^2*zeta)/8 - eta^2/8 - zeta^2/8 - (eta*ksi*zeta)/4 + 1/8, ksi/4 - (eta*ksi)/4 - (eta*zeta)/8 + (ksi*zeta)/4 - (eta*zeta^2)/8 + (eta^2*zeta)/8 + eta^2/8 + zeta^2/8 - (eta*ksi*zeta)/4 - 1/8, ksi/4 + (eta*ksi)/4 + (eta*zeta)/8 + (ksi*zeta)/4 + (eta*zeta^2)/8 + (eta^2*zeta)/8 + eta^2/8 + zeta^2/8 + (eta*ksi*zeta)/4 - 1/8, ksi/4 + (eta*ksi)/4 - (eta*zeta)/8 + (ksi*zeta)/4 - (eta*zeta^2)/8 - (eta^2*zeta)/8 - eta^2/8 - zeta^2/8 + (eta*ksi*zeta)/4 + 1/8, (eta*ksi)/2 - ksi/2 + (ksi*zeta)/2 - (eta*ksi*zeta)/2,               (eta^2*zeta)/4 - zeta/4 - eta^2/4 + 1/4, (ksi*zeta)/2 - (eta*ksi)/2 - ksi/2 + (eta*ksi*zeta)/2,               zeta/4 - (eta^2*zeta)/4 + eta^2/4 - 1/4,                 eta/4 - (eta*zeta^2)/4 + zeta^2/4 - 1/4,                 (eta*zeta^2)/4 - eta/4 - zeta^2/4 + 1/4,                   eta/4 - (eta*zeta^2)/4 - zeta^2/4 + 1/4,                 (eta*zeta^2)/4 - eta/4 + zeta^2/4 - 1/4, (eta*ksi)/2 - ksi/2 - (ksi*zeta)/2 + (eta*ksi*zeta)/2,                 zeta/4 - (eta^2*zeta)/4 - eta^2/4 + 1/4, - ksi/2 - (eta*ksi)/2 - (ksi*zeta)/2 - (eta*ksi*zeta)/2,               (eta^2*zeta)/4 - zeta/4 + eta^2/4 - 1/4
            eta/4 - (eta*ksi)/4 - (eta*zeta)/4 - (ksi*zeta)/8 + (ksi*zeta^2)/8 + (ksi^2*zeta)/8 - ksi^2/8 - zeta^2/8 + (eta*ksi*zeta)/4 + 1/8, eta/4 + (eta*ksi)/4 - (eta*zeta)/4 + (ksi*zeta)/8 - (ksi*zeta^2)/8 + (ksi^2*zeta)/8 - ksi^2/8 - zeta^2/8 - (eta*ksi*zeta)/4 + 1/8, eta/4 + (eta*ksi)/4 - (eta*zeta)/4 - (ksi*zeta)/8 + (ksi*zeta^2)/8 - (ksi^2*zeta)/8 + ksi^2/8 + zeta^2/8 - (eta*ksi*zeta)/4 - 1/8, eta/4 - (eta*ksi)/4 - (eta*zeta)/4 + (ksi*zeta)/8 - (ksi*zeta^2)/8 - (ksi^2*zeta)/8 + ksi^2/8 + zeta^2/8 + (eta*ksi*zeta)/4 - 1/8, eta/4 - (eta*ksi)/4 + (eta*zeta)/4 + (ksi*zeta)/8 + (ksi*zeta^2)/8 - (ksi^2*zeta)/8 - ksi^2/8 - zeta^2/8 - (eta*ksi*zeta)/4 + 1/8, eta/4 + (eta*ksi)/4 + (eta*zeta)/4 - (ksi*zeta)/8 - (ksi*zeta^2)/8 - (ksi^2*zeta)/8 - ksi^2/8 - zeta^2/8 + (eta*ksi*zeta)/4 + 1/8, eta/4 + (eta*ksi)/4 + (eta*zeta)/4 + (ksi*zeta)/8 + (ksi*zeta^2)/8 + (ksi^2*zeta)/8 + ksi^2/8 + zeta^2/8 + (eta*ksi*zeta)/4 - 1/8, eta/4 - (eta*ksi)/4 + (eta*zeta)/4 - (ksi*zeta)/8 - (ksi*zeta^2)/8 + (ksi^2*zeta)/8 + ksi^2/8 + zeta^2/8 - (eta*ksi*zeta)/4 - 1/8,               zeta/4 - (ksi^2*zeta)/4 + ksi^2/4 - 1/4, (eta*zeta)/2 - (eta*ksi)/2 - eta/2 + (eta*ksi*zeta)/2,               (ksi^2*zeta)/4 - zeta/4 - ksi^2/4 + 1/4, (eta*ksi)/2 - eta/2 + (eta*zeta)/2 - (eta*ksi*zeta)/2,                 ksi/4 - (ksi*zeta^2)/4 + zeta^2/4 - 1/4,                 (ksi*zeta^2)/4 - ksi/4 + zeta^2/4 - 1/4,                   ksi/4 - (ksi*zeta^2)/4 - zeta^2/4 + 1/4,                 (ksi*zeta^2)/4 - ksi/4 - zeta^2/4 + 1/4,               (ksi^2*zeta)/4 - zeta/4 + ksi^2/4 - 1/4, - eta/2 - (eta*ksi)/2 - (eta*zeta)/2 - (eta*ksi*zeta)/2,                 zeta/4 - (ksi^2*zeta)/4 - ksi^2/4 + 1/4, (eta*ksi)/2 - eta/2 - (eta*zeta)/2 + (eta*ksi*zeta)/2
            zeta/4 - (eta*ksi)/8 - (eta*zeta)/4 - (ksi*zeta)/4 + (eta*ksi^2)/8 + (eta^2*ksi)/8 - eta^2/8 - ksi^2/8 + (eta*ksi*zeta)/4 + 1/8,   zeta/4 + (eta*ksi)/8 - (eta*zeta)/4 + (ksi*zeta)/4 + (eta*ksi^2)/8 - (eta^2*ksi)/8 - eta^2/8 - ksi^2/8 - (eta*ksi*zeta)/4 + 1/8,   zeta/4 - (eta*ksi)/8 + (eta*zeta)/4 + (ksi*zeta)/4 - (eta*ksi^2)/8 - (eta^2*ksi)/8 - eta^2/8 - ksi^2/8 + (eta*ksi*zeta)/4 + 1/8,   zeta/4 + (eta*ksi)/8 + (eta*zeta)/4 - (ksi*zeta)/4 - (eta*ksi^2)/8 + (eta^2*ksi)/8 - eta^2/8 - ksi^2/8 - (eta*ksi*zeta)/4 + 1/8,   zeta/4 + (eta*ksi)/8 - (eta*zeta)/4 - (ksi*zeta)/4 - (eta*ksi^2)/8 - (eta^2*ksi)/8 + eta^2/8 + ksi^2/8 + (eta*ksi*zeta)/4 - 1/8,   zeta/4 - (eta*ksi)/8 - (eta*zeta)/4 + (ksi*zeta)/4 - (eta*ksi^2)/8 + (eta^2*ksi)/8 + eta^2/8 + ksi^2/8 - (eta*ksi*zeta)/4 - 1/8,   zeta/4 + (eta*ksi)/8 + (eta*zeta)/4 + (ksi*zeta)/4 + (eta*ksi^2)/8 + (eta^2*ksi)/8 + eta^2/8 + ksi^2/8 + (eta*ksi*zeta)/4 - 1/8,   zeta/4 - (eta*ksi)/8 + (eta*zeta)/4 - (ksi*zeta)/4 + (eta*ksi^2)/8 - (eta^2*ksi)/8 + eta^2/8 + ksi^2/8 - (eta*ksi*zeta)/4 - 1/8,                 eta/4 - (eta*ksi^2)/4 + ksi^2/4 - 1/4,                 (eta^2*ksi)/4 - ksi/4 + eta^2/4 - 1/4,                 (eta*ksi^2)/4 - eta/4 + ksi^2/4 - 1/4,                 ksi/4 - (eta^2*ksi)/4 + eta^2/4 - 1/4, (eta*zeta)/2 - zeta/2 + (ksi*zeta)/2 - (eta*ksi*zeta)/2, (eta*zeta)/2 - zeta/2 - (ksi*zeta)/2 + (eta*ksi*zeta)/2, - zeta/2 - (eta*zeta)/2 - (ksi*zeta)/2 - (eta*ksi*zeta)/2, (ksi*zeta)/2 - (eta*zeta)/2 - zeta/2 + (eta*ksi*zeta)/2,                 (eta*ksi^2)/4 - eta/4 - ksi^2/4 + 1/4,                   ksi/4 - (eta^2*ksi)/4 - eta^2/4 + 1/4,                   eta/4 - (eta*ksi^2)/4 - ksi^2/4 + 1/4,                 (eta^2*ksi)/4 - ksi/4 - eta^2/4 + 1/4];
        
        % Derivadas de x,y,z respecto de ksi, eta, zeta
        jac = dN*nodesEle;
        
        dNxyz = jac\dN;
        
        B = zeros(size(C,2),nDofNod*nNodEle);
        B(1,1:3:(length(B)-2)) = dNxyz(1,:);
        B(2,2:3:(length(B)-1)) = dNxyz(2,:);
        B(3,3:3:length(B))     = dNxyz(3,:);
        B(4,1:3:(length(B)-2)) = dNxyz(2,:);
        B(4,2:3:(length(B)-1)) = dNxyz(1,:);
        B(5,2:3:(length(B)-1)) = dNxyz(3,:);
        B(5,3:3:length(B))     = dNxyz(2,:);
        B(6,1:3:(length(B)-2)) = dNxyz(3,:);
        B(6,3:3:length(B))     = dNxyz(1,:);
        
        Ke = Ke + (B'*C*B*wpg(ipg)*det(jac));
        Fe = Fe + (B'*C*epsilon*wpg(ipg)*det(jac));
    end
    eleDofs = nodeDofs(elements(iele,:),:);
    eleDofs = reshape(eleDofs',[],1);
    K(eleDofs,eleDofs) = K(eleDofs,eleDofs) + Ke;
    F(eleDofs) = F(eleDofs) + Fe;
end



%% Calculo de cargas consistentes
switch Caso
    case 'Normal'
        if T==0
        P=-10;
        Area=LadoH*Espesor;
        AreaElemento=Area/(ElementosEspesor*ElementosH);
        p=P/Area*AreaElemento;
        for iele = 1:(ElementosEspesor*ElementosH)           
            R(elements(ElementosLargo*iele,[2 3 6 7]),3)=R(elements(ElementosLargo*iele,[2 3 6 7]),3)-1/12*p;
            R(elements(ElementosLargo*iele,[10 14 15 18]),3)=R(elements(ElementosLargo*iele,[10 14 15 18]),3)+1/3*p;
        end
        end
end







%% Solver
% Reduccion Matriz
isFixed = logical(reshape(bc',[],1));
isFree = ~isFixed;

Rr = reshape(R',[],1);
K=sparse(K);
Rr=sparse(Rr)+F;
% Solver
Dr = K(isFree,isFree)\Rr(isFree);

% Reconstrucciï¿½n
D = zeros(nDofTot,1);
D(isFree) = D(isFree) + Dr;
Dnodes=nodes+[D(1:3:length(D)) D(2:3:length(D)) D(3:3:length(D))];
% scatter3(nodes(:,1),nodes(:,2),nodes(:,3))
% hold on
% scatter3(Dnodes(:,1),Dnodes(:,2),Dnodes(:,3),'r')
% daspect([1 1 1]);
% figure
% meshplot(elements,nodes,'b')
% hold on
% meshplot(elements,Dnodes,'r')
% Reacciones
Rv = K(isFixed,isFree)*D(isFree);
reacciones = nan(nDofTot,1);
reacciones(isFixed) = Rv;
reacciones = (reshape(reacciones,nDofNod,[]))';


%% Tensiones en Puntos de Gauss Full 27
StressGauss27 = zeros(nel,npg,6);
% xgauss=[ones(4,1);-ones(4,1)];
%  ygauss=repmat(xgauss([1 2 5 6]),2,1);
%  zgauss=repmat([-1 1]',4,1);
%  upg=a*[xgauss ygauss zgauss]; % Ubicaciones puntos de Gauss
%  wpg=ones(8,1);
%  npg = size(upg,1); % Numero de puntos de Gauss
for iele = 1:nel
    nodesEle = nodes(elements(iele,:),:);
    if max(nodes(elements(iele,:),3))<= Espesor1
        C=C1;
    else
        C=C2;
    end
    for ipg = 1:npg
        % Punto de Gauss
        ksi = upg(ipg,1);
        eta = upg(ipg,2);
        zeta = upg(ipg,3);
        % Derivadas de las funciones de forma respecto de ksi, eta
        %         [N dN]=FuncionFormaH20(ksi,eta,zeta);
        dN=[ ksi/4 - (eta*ksi)/4 - (eta*zeta)/8 - (ksi*zeta)/4 + (eta*zeta^2)/8 + (eta^2*zeta)/8 - eta^2/8 - zeta^2/8 + (eta*ksi*zeta)/4 + 1/8, ksi/4 - (eta*ksi)/4 + (eta*zeta)/8 - (ksi*zeta)/4 - (eta*zeta^2)/8 - (eta^2*zeta)/8 + eta^2/8 + zeta^2/8 + (eta*ksi*zeta)/4 - 1/8, ksi/4 + (eta*ksi)/4 - (eta*zeta)/8 - (ksi*zeta)/4 + (eta*zeta^2)/8 - (eta^2*zeta)/8 + eta^2/8 + zeta^2/8 - (eta*ksi*zeta)/4 - 1/8, ksi/4 + (eta*ksi)/4 + (eta*zeta)/8 - (ksi*zeta)/4 - (eta*zeta^2)/8 + (eta^2*zeta)/8 - eta^2/8 - zeta^2/8 - (eta*ksi*zeta)/4 + 1/8, ksi/4 - (eta*ksi)/4 + (eta*zeta)/8 + (ksi*zeta)/4 + (eta*zeta^2)/8 - (eta^2*zeta)/8 - eta^2/8 - zeta^2/8 - (eta*ksi*zeta)/4 + 1/8, ksi/4 - (eta*ksi)/4 - (eta*zeta)/8 + (ksi*zeta)/4 - (eta*zeta^2)/8 + (eta^2*zeta)/8 + eta^2/8 + zeta^2/8 - (eta*ksi*zeta)/4 - 1/8, ksi/4 + (eta*ksi)/4 + (eta*zeta)/8 + (ksi*zeta)/4 + (eta*zeta^2)/8 + (eta^2*zeta)/8 + eta^2/8 + zeta^2/8 + (eta*ksi*zeta)/4 - 1/8, ksi/4 + (eta*ksi)/4 - (eta*zeta)/8 + (ksi*zeta)/4 - (eta*zeta^2)/8 - (eta^2*zeta)/8 - eta^2/8 - zeta^2/8 + (eta*ksi*zeta)/4 + 1/8, (eta*ksi)/2 - ksi/2 + (ksi*zeta)/2 - (eta*ksi*zeta)/2,               (eta^2*zeta)/4 - zeta/4 - eta^2/4 + 1/4, (ksi*zeta)/2 - (eta*ksi)/2 - ksi/2 + (eta*ksi*zeta)/2,               zeta/4 - (eta^2*zeta)/4 + eta^2/4 - 1/4,                 eta/4 - (eta*zeta^2)/4 + zeta^2/4 - 1/4,                 (eta*zeta^2)/4 - eta/4 - zeta^2/4 + 1/4,                   eta/4 - (eta*zeta^2)/4 - zeta^2/4 + 1/4,                 (eta*zeta^2)/4 - eta/4 + zeta^2/4 - 1/4, (eta*ksi)/2 - ksi/2 - (ksi*zeta)/2 + (eta*ksi*zeta)/2,                 zeta/4 - (eta^2*zeta)/4 - eta^2/4 + 1/4, - ksi/2 - (eta*ksi)/2 - (ksi*zeta)/2 - (eta*ksi*zeta)/2,               (eta^2*zeta)/4 - zeta/4 + eta^2/4 - 1/4
            eta/4 - (eta*ksi)/4 - (eta*zeta)/4 - (ksi*zeta)/8 + (ksi*zeta^2)/8 + (ksi^2*zeta)/8 - ksi^2/8 - zeta^2/8 + (eta*ksi*zeta)/4 + 1/8, eta/4 + (eta*ksi)/4 - (eta*zeta)/4 + (ksi*zeta)/8 - (ksi*zeta^2)/8 + (ksi^2*zeta)/8 - ksi^2/8 - zeta^2/8 - (eta*ksi*zeta)/4 + 1/8, eta/4 + (eta*ksi)/4 - (eta*zeta)/4 - (ksi*zeta)/8 + (ksi*zeta^2)/8 - (ksi^2*zeta)/8 + ksi^2/8 + zeta^2/8 - (eta*ksi*zeta)/4 - 1/8, eta/4 - (eta*ksi)/4 - (eta*zeta)/4 + (ksi*zeta)/8 - (ksi*zeta^2)/8 - (ksi^2*zeta)/8 + ksi^2/8 + zeta^2/8 + (eta*ksi*zeta)/4 - 1/8, eta/4 - (eta*ksi)/4 + (eta*zeta)/4 + (ksi*zeta)/8 + (ksi*zeta^2)/8 - (ksi^2*zeta)/8 - ksi^2/8 - zeta^2/8 - (eta*ksi*zeta)/4 + 1/8, eta/4 + (eta*ksi)/4 + (eta*zeta)/4 - (ksi*zeta)/8 - (ksi*zeta^2)/8 - (ksi^2*zeta)/8 - ksi^2/8 - zeta^2/8 + (eta*ksi*zeta)/4 + 1/8, eta/4 + (eta*ksi)/4 + (eta*zeta)/4 + (ksi*zeta)/8 + (ksi*zeta^2)/8 + (ksi^2*zeta)/8 + ksi^2/8 + zeta^2/8 + (eta*ksi*zeta)/4 - 1/8, eta/4 - (eta*ksi)/4 + (eta*zeta)/4 - (ksi*zeta)/8 - (ksi*zeta^2)/8 + (ksi^2*zeta)/8 + ksi^2/8 + zeta^2/8 - (eta*ksi*zeta)/4 - 1/8,               zeta/4 - (ksi^2*zeta)/4 + ksi^2/4 - 1/4, (eta*zeta)/2 - (eta*ksi)/2 - eta/2 + (eta*ksi*zeta)/2,               (ksi^2*zeta)/4 - zeta/4 - ksi^2/4 + 1/4, (eta*ksi)/2 - eta/2 + (eta*zeta)/2 - (eta*ksi*zeta)/2,                 ksi/4 - (ksi*zeta^2)/4 + zeta^2/4 - 1/4,                 (ksi*zeta^2)/4 - ksi/4 + zeta^2/4 - 1/4,                   ksi/4 - (ksi*zeta^2)/4 - zeta^2/4 + 1/4,                 (ksi*zeta^2)/4 - ksi/4 - zeta^2/4 + 1/4,               (ksi^2*zeta)/4 - zeta/4 + ksi^2/4 - 1/4, - eta/2 - (eta*ksi)/2 - (eta*zeta)/2 - (eta*ksi*zeta)/2,                 zeta/4 - (ksi^2*zeta)/4 - ksi^2/4 + 1/4, (eta*ksi)/2 - eta/2 - (eta*zeta)/2 + (eta*ksi*zeta)/2
            zeta/4 - (eta*ksi)/8 - (eta*zeta)/4 - (ksi*zeta)/4 + (eta*ksi^2)/8 + (eta^2*ksi)/8 - eta^2/8 - ksi^2/8 + (eta*ksi*zeta)/4 + 1/8,   zeta/4 + (eta*ksi)/8 - (eta*zeta)/4 + (ksi*zeta)/4 + (eta*ksi^2)/8 - (eta^2*ksi)/8 - eta^2/8 - ksi^2/8 - (eta*ksi*zeta)/4 + 1/8,   zeta/4 - (eta*ksi)/8 + (eta*zeta)/4 + (ksi*zeta)/4 - (eta*ksi^2)/8 - (eta^2*ksi)/8 - eta^2/8 - ksi^2/8 + (eta*ksi*zeta)/4 + 1/8,   zeta/4 + (eta*ksi)/8 + (eta*zeta)/4 - (ksi*zeta)/4 - (eta*ksi^2)/8 + (eta^2*ksi)/8 - eta^2/8 - ksi^2/8 - (eta*ksi*zeta)/4 + 1/8,   zeta/4 + (eta*ksi)/8 - (eta*zeta)/4 - (ksi*zeta)/4 - (eta*ksi^2)/8 - (eta^2*ksi)/8 + eta^2/8 + ksi^2/8 + (eta*ksi*zeta)/4 - 1/8,   zeta/4 - (eta*ksi)/8 - (eta*zeta)/4 + (ksi*zeta)/4 - (eta*ksi^2)/8 + (eta^2*ksi)/8 + eta^2/8 + ksi^2/8 - (eta*ksi*zeta)/4 - 1/8,   zeta/4 + (eta*ksi)/8 + (eta*zeta)/4 + (ksi*zeta)/4 + (eta*ksi^2)/8 + (eta^2*ksi)/8 + eta^2/8 + ksi^2/8 + (eta*ksi*zeta)/4 - 1/8,   zeta/4 - (eta*ksi)/8 + (eta*zeta)/4 - (ksi*zeta)/4 + (eta*ksi^2)/8 - (eta^2*ksi)/8 + eta^2/8 + ksi^2/8 - (eta*ksi*zeta)/4 - 1/8,                 eta/4 - (eta*ksi^2)/4 + ksi^2/4 - 1/4,                 (eta^2*ksi)/4 - ksi/4 + eta^2/4 - 1/4,                 (eta*ksi^2)/4 - eta/4 + ksi^2/4 - 1/4,                 ksi/4 - (eta^2*ksi)/4 + eta^2/4 - 1/4, (eta*zeta)/2 - zeta/2 + (ksi*zeta)/2 - (eta*ksi*zeta)/2, (eta*zeta)/2 - zeta/2 - (ksi*zeta)/2 + (eta*ksi*zeta)/2, - zeta/2 - (eta*zeta)/2 - (ksi*zeta)/2 - (eta*ksi*zeta)/2, (ksi*zeta)/2 - (eta*zeta)/2 - zeta/2 + (eta*ksi*zeta)/2,                 (eta*ksi^2)/4 - eta/4 - ksi^2/4 + 1/4,                   ksi/4 - (eta^2*ksi)/4 - eta^2/4 + 1/4,                   eta/4 - (eta*ksi^2)/4 - ksi^2/4 + 1/4,                 (eta^2*ksi)/4 - ksi/4 - eta^2/4 + 1/4];
        
        % Derivadas de x,y, respecto de ksi, eta
        jac = dN*nodesEle;
        % Derivadas de las funciones de forma respecto de x,y.
        dNxyz = jac\dN;          % dNxy = inv(jac)*dN
        
        B = zeros(size(C,2),nDofNod*nNodEle);
        B(1,1:3:(length(B)-2)) = dNxyz(1,:);
        B(2,2:3:(length(B)-1)) = dNxyz(2,:);
        B(3,3:3:length(B))     = dNxyz(3,:);
        B(4,1:3:(length(B)-2)) = dNxyz(2,:);
        B(4,2:3:(length(B)-1)) = dNxyz(1,:);
        B(5,2:3:(length(B)-1)) = dNxyz(3,:);
        B(5,3:3:length(B))     = dNxyz(2,:);
        B(6,1:3:(length(B)-2)) = dNxyz(3,:);
        B(6,3:3:length(B))     = dNxyz(1,:);
        
        eleDofs = nodeDofs(elements(iele,:),:);
        eleDofs = reshape(eleDofs',[],1);
        StressGauss27(iele,ipg,:) = C*B*D(eleDofs)-C*epsilon;
    end
end


%% Extrapolo a Nodos las tensiones en puntos de gauss Full
Nodos=[-1    -1   -1
        1    -1   -1
        1     1   -1
       -1     1   -1
       -1    -1    1
        1    -1    1
        1     1    1
       -1     1    1
        0    -1   -1
        1     0   -1
        0     1   -1
       -1     0   -1
       -1    -1    0
        1    -1    0
        1     1    0
       -1     1    0
        0    -1    1
        1     0    1
        0     1    1
       -1     0    1]*(5/3)^0.5;
Func= @(ksi,eta,zeta) [ (eta.^2.*ksi.^2.*zeta.^2)/8 - (eta.^2.*ksi.^2.*zeta)/8 - (eta.^2.*ksi.*zeta.^2)/8 + (eta.^2.*ksi.*zeta)/8 - (eta.*ksi.^2.*zeta.^2)/8 + (eta.*ksi.^2.*zeta)/8 + (eta.*ksi.*zeta.^2)/8 - (eta.*ksi.*zeta)/8, (eta.^2.*ksi.^2.*zeta.^2)/8 - (eta.^2.*ksi.^2.*zeta)/8 + (eta.^2.*ksi.*zeta.^2)/8 - (eta.^2.*ksi.*zeta)/8 - (eta.*ksi.^2.*zeta.^2)/8 + (eta.*ksi.^2.*zeta)/8 - (eta.*ksi.*zeta.^2)/8 + (eta.*ksi.*zeta)/8, (eta.^2.*ksi.^2.*zeta.^2)/8 - (eta.^2.*ksi.^2.*zeta)/8 + (eta.^2.*ksi.*zeta.^2)/8 - (eta.^2.*ksi.*zeta)/8 + (eta.*ksi.^2.*zeta.^2)/8 - (eta.*ksi.^2.*zeta)/8 + (eta.*ksi.*zeta.^2)/8 - (eta.*ksi.*zeta)/8, (eta.^2.*ksi.^2.*zeta.^2)/8 - (eta.^2.*ksi.^2.*zeta)/8 - (eta.^2.*ksi.*zeta.^2)/8 + (eta.^2.*ksi.*zeta)/8 + (eta.*ksi.^2.*zeta.^2)/8 - (eta.*ksi.^2.*zeta)/8 - (eta.*ksi.*zeta.^2)/8 + (eta.*ksi.*zeta)/8, (eta.^2.*ksi.^2.*zeta.^2)/8 + (eta.^2.*ksi.^2.*zeta)/8 - (eta.^2.*ksi.*zeta.^2)/8 - (eta.^2.*ksi.*zeta)/8 - (eta.*ksi.^2.*zeta.^2)/8 - (eta.*ksi.^2.*zeta)/8 + (eta.*ksi.*zeta.^2)/8 + (eta.*ksi.*zeta)/8, (eta.^2.*ksi.^2.*zeta.^2)/8 + (eta.^2.*ksi.^2.*zeta)/8 + (eta.^2.*ksi.*zeta.^2)/8 + (eta.^2.*ksi.*zeta)/8 - (eta.*ksi.^2.*zeta.^2)/8 - (eta.*ksi.^2.*zeta)/8 - (eta.*ksi.*zeta.^2)/8 - (eta.*ksi.*zeta)/8, (eta.^2.*ksi.^2.*zeta.^2)/8 + (eta.^2.*ksi.^2.*zeta)/8 + (eta.^2.*ksi.*zeta.^2)/8 + (eta.^2.*ksi.*zeta)/8 + (eta.*ksi.^2.*zeta.^2)/8 + (eta.*ksi.^2.*zeta)/8 + (eta.*ksi.*zeta.^2)/8 + (eta.*ksi.*zeta)/8, (eta.^2.*ksi.^2.*zeta.^2)/8 + (eta.^2.*ksi.^2.*zeta)/8 - (eta.^2.*ksi.*zeta.^2)/8 - (eta.^2.*ksi.*zeta)/8 + (eta.*ksi.^2.*zeta.^2)/8 + (eta.*ksi.^2.*zeta)/8 - (eta.*ksi.*zeta.^2)/8 - (eta.*ksi.*zeta)/8, - (eta.^2.*ksi.^2.*zeta.^2)/4 + (eta.^2.*ksi.^2.*zeta)/4 + (eta.^2.*zeta.^2)/4 - (eta.^2.*zeta)/4 + (eta.*ksi.^2.*zeta.^2)/4 - (eta.*ksi.^2.*zeta)/4 - (eta.*zeta.^2)/4 + (eta.*zeta)/4, - (eta.^2.*ksi.^2.*zeta.^2)/4 + (eta.^2.*ksi.^2.*zeta)/4 - (eta.^2.*ksi.*zeta.^2)/4 + (eta.^2.*ksi.*zeta)/4 + (ksi.^2.*zeta.^2)/4 - (ksi.^2.*zeta)/4 + (ksi.*zeta.^2)/4 - (ksi.*zeta)/4, - (eta.^2.*ksi.^2.*zeta.^2)/4 + (eta.^2.*ksi.^2.*zeta)/4 + (eta.^2.*zeta.^2)/4 - (eta.^2.*zeta)/4 - (eta.*ksi.^2.*zeta.^2)/4 + (eta.*ksi.^2.*zeta)/4 + (eta.*zeta.^2)/4 - (eta.*zeta)/4, - (eta.^2.*ksi.^2.*zeta.^2)/4 + (eta.^2.*ksi.^2.*zeta)/4 + (eta.^2.*ksi.*zeta.^2)/4 - (eta.^2.*ksi.*zeta)/4 + (ksi.^2.*zeta.^2)/4 - (ksi.^2.*zeta)/4 - (ksi.*zeta.^2)/4 + (ksi.*zeta)/4, - (eta.^2.*ksi.^2.*zeta.^2)/4 + (eta.^2.*ksi.^2)/4 + (eta.^2.*ksi.*zeta.^2)/4 - (eta.^2.*ksi)/4 + (eta.*ksi.^2.*zeta.^2)/4 - (eta.*ksi.^2)/4 - (eta.*ksi.*zeta.^2)/4 + (eta.*ksi)/4, - (eta.^2.*ksi.^2.*zeta.^2)/4 + (eta.^2.*ksi.^2)/4 - (eta.^2.*ksi.*zeta.^2)/4 + (eta.^2.*ksi)/4 + (eta.*ksi.^2.*zeta.^2)/4 - (eta.*ksi.^2)/4 + (eta.*ksi.*zeta.^2)/4 - (eta.*ksi)/4, (eta.*ksi)/4 + (eta.*ksi.^2)/4 + (eta.^2.*ksi)/4 + (eta.^2.*ksi.^2)/4 - (eta.*ksi.*zeta.^2)/4 - (eta.*ksi.^2.*zeta.^2)/4 - (eta.^2.*ksi.*zeta.^2)/4 - (eta.^2.*ksi.^2.*zeta.^2)/4, - (eta.^2.*ksi.^2.*zeta.^2)/4 + (eta.^2.*ksi.^2)/4 + (eta.^2.*ksi.*zeta.^2)/4 - (eta.^2.*ksi)/4 - (eta.*ksi.^2.*zeta.^2)/4 + (eta.*ksi.^2)/4 + (eta.*ksi.*zeta.^2)/4 - (eta.*ksi)/4, - (eta.^2.*ksi.^2.*zeta.^2)/4 - (eta.^2.*ksi.^2.*zeta)/4 + (eta.^2.*zeta.^2)/4 + (eta.^2.*zeta)/4 + (eta.*ksi.^2.*zeta.^2)/4 + (eta.*ksi.^2.*zeta)/4 - (eta.*zeta.^2)/4 - (eta.*zeta)/4, (ksi.^2.*zeta.^2)/4 + (ksi.*zeta)/4 + (ksi.*zeta.^2)/4 + (ksi.^2.*zeta)/4 - (eta.^2.*ksi.*zeta)/4 - (eta.^2.*ksi.*zeta.^2)/4 - (eta.^2.*ksi.^2.*zeta)/4 - (eta.^2.*ksi.^2.*zeta.^2)/4, (eta.*zeta)/4 + (eta.*zeta.^2)/4 + (eta.^2.*zeta)/4 + (eta.^2.*zeta.^2)/4 - (eta.*ksi.^2.*zeta)/4 - (eta.*ksi.^2.*zeta.^2)/4 - (eta.^2.*ksi.^2.*zeta)/4 - (eta.^2.*ksi.^2.*zeta.^2)/4, - (eta.^2.*ksi.^2.*zeta.^2)/4 - (eta.^2.*ksi.^2.*zeta)/4 + (eta.^2.*ksi.*zeta.^2)/4 + (eta.^2.*ksi.*zeta)/4 + (ksi.^2.*zeta.^2)/4 + (ksi.^2.*zeta)/4 - (ksi.*zeta.^2)/4 - (ksi.*zeta)/4, (eta.^2.*ksi.^2.*zeta.^2)/2 - (eta.^2.*ksi.^2.*zeta)/2 - (eta.^2.*zeta.^2)/2 + (eta.^2.*zeta)/2 - (ksi.^2.*zeta.^2)/2 + (ksi.^2.*zeta)/2 + zeta.^2/2 - zeta/2, (eta.^2.*ksi.^2.*zeta.^2)/2 - (eta.^2.*ksi.^2)/2 - (eta.^2.*ksi.*zeta.^2)/2 + (eta.^2.*ksi)/2 - (ksi.^2.*zeta.^2)/2 + ksi.^2/2 + (ksi.*zeta.^2)/2 - ksi/2, (eta.^2.*ksi.^2.*zeta.^2)/2 - (eta.^2.*ksi.^2)/2 - (eta.^2.*zeta.^2)/2 + eta.^2/2 - (eta.*ksi.^2.*zeta.^2)/2 + (eta.*ksi.^2)/2 + (eta.*zeta.^2)/2 - eta/2, (eta.^2.*ksi.^2.*zeta.^2)/2 - (eta.^2.*ksi.^2)/2 + (eta.^2.*ksi.*zeta.^2)/2 - (eta.^2.*ksi)/2 - (ksi.^2.*zeta.^2)/2 + ksi.^2/2 - (ksi.*zeta.^2)/2 + ksi/2, (eta.^2.*ksi.^2.*zeta.^2)/2 - (eta.^2.*ksi.^2)/2 - (eta.^2.*zeta.^2)/2 + eta.^2/2 + (eta.*ksi.^2.*zeta.^2)/2 - (eta.*ksi.^2)/2 - (eta.*zeta.^2)/2 + eta/2, (eta.^2.*ksi.^2.*zeta.^2)/2 + (eta.^2.*ksi.^2.*zeta)/2 - (eta.^2.*zeta.^2)/2 - (eta.^2.*zeta)/2 - (ksi.^2.*zeta.^2)/2 - (ksi.^2.*zeta)/2 + zeta.^2/2 + zeta/2, - eta.^2.*ksi.^2.*zeta.^2 + eta.^2.*ksi.^2 + eta.^2.*zeta.^2 - eta.^2 + ksi.^2.*zeta.^2 - ksi.^2 - zeta.^2 + 1];
N= Func(Nodos(:,1),Nodos(:,2),Nodos(:,3));
stressExt27= zeros(nel,nNodEle,6);
AuxStressGauss27=permute(StressGauss27,[2 3 1]);
for iele = 1:nel
    stressExt27(iele,:,:) = N*AuxStressGauss27(:,:,iele);
end

%% Tensiones en Puntos de Gauss Reducido 8

xyzgauss=[-1    -1   -1
        1    -1   -1
        1     1   -1
       -1     1   -1
       -1    -1    1
        1    -1    1
        1     1    1
       -1     1    1];
upg=[xyzgauss]/3^0.5; % Ubicaciones puntos de Gauss
wpg=ones(8,1);
npg = size(upg,1); % Numero de puntos de Gauss
StressGauss8 = zeros(nel,npg,6);
for iele = 1:nel
    nodesEle = nodes(elements(iele,:),:);
    if max(nodes(elements(iele,:),3))<= Espesor1
        C=C1;
        epsilon=alfa1*T*[1;1;1;0;0;0];
    else
        C=C2;
        epsilon=alfa2*T*[1;1;1;0;0;0];
    end
    for ipg = 1:npg
        % Punto de Gauss
        ksi = upg(ipg,1);
        eta = upg(ipg,2);
        zeta = upg(ipg,3);
        % Derivadas de las funciones de forma respecto de ksi, eta
        %         [N dN]=FuncionFormaH20(ksi,eta,zeta);
        dN=[ ksi/4 - (eta*ksi)/4 - (eta*zeta)/8 - (ksi*zeta)/4 + (eta*zeta^2)/8 + (eta^2*zeta)/8 - eta^2/8 - zeta^2/8 + (eta*ksi*zeta)/4 + 1/8, ksi/4 - (eta*ksi)/4 + (eta*zeta)/8 - (ksi*zeta)/4 - (eta*zeta^2)/8 - (eta^2*zeta)/8 + eta^2/8 + zeta^2/8 + (eta*ksi*zeta)/4 - 1/8, ksi/4 + (eta*ksi)/4 - (eta*zeta)/8 - (ksi*zeta)/4 + (eta*zeta^2)/8 - (eta^2*zeta)/8 + eta^2/8 + zeta^2/8 - (eta*ksi*zeta)/4 - 1/8, ksi/4 + (eta*ksi)/4 + (eta*zeta)/8 - (ksi*zeta)/4 - (eta*zeta^2)/8 + (eta^2*zeta)/8 - eta^2/8 - zeta^2/8 - (eta*ksi*zeta)/4 + 1/8, ksi/4 - (eta*ksi)/4 + (eta*zeta)/8 + (ksi*zeta)/4 + (eta*zeta^2)/8 - (eta^2*zeta)/8 - eta^2/8 - zeta^2/8 - (eta*ksi*zeta)/4 + 1/8, ksi/4 - (eta*ksi)/4 - (eta*zeta)/8 + (ksi*zeta)/4 - (eta*zeta^2)/8 + (eta^2*zeta)/8 + eta^2/8 + zeta^2/8 - (eta*ksi*zeta)/4 - 1/8, ksi/4 + (eta*ksi)/4 + (eta*zeta)/8 + (ksi*zeta)/4 + (eta*zeta^2)/8 + (eta^2*zeta)/8 + eta^2/8 + zeta^2/8 + (eta*ksi*zeta)/4 - 1/8, ksi/4 + (eta*ksi)/4 - (eta*zeta)/8 + (ksi*zeta)/4 - (eta*zeta^2)/8 - (eta^2*zeta)/8 - eta^2/8 - zeta^2/8 + (eta*ksi*zeta)/4 + 1/8, (eta*ksi)/2 - ksi/2 + (ksi*zeta)/2 - (eta*ksi*zeta)/2,               (eta^2*zeta)/4 - zeta/4 - eta^2/4 + 1/4, (ksi*zeta)/2 - (eta*ksi)/2 - ksi/2 + (eta*ksi*zeta)/2,               zeta/4 - (eta^2*zeta)/4 + eta^2/4 - 1/4,                 eta/4 - (eta*zeta^2)/4 + zeta^2/4 - 1/4,                 (eta*zeta^2)/4 - eta/4 - zeta^2/4 + 1/4,                   eta/4 - (eta*zeta^2)/4 - zeta^2/4 + 1/4,                 (eta*zeta^2)/4 - eta/4 + zeta^2/4 - 1/4, (eta*ksi)/2 - ksi/2 - (ksi*zeta)/2 + (eta*ksi*zeta)/2,                 zeta/4 - (eta^2*zeta)/4 - eta^2/4 + 1/4, - ksi/2 - (eta*ksi)/2 - (ksi*zeta)/2 - (eta*ksi*zeta)/2,               (eta^2*zeta)/4 - zeta/4 + eta^2/4 - 1/4
            eta/4 - (eta*ksi)/4 - (eta*zeta)/4 - (ksi*zeta)/8 + (ksi*zeta^2)/8 + (ksi^2*zeta)/8 - ksi^2/8 - zeta^2/8 + (eta*ksi*zeta)/4 + 1/8, eta/4 + (eta*ksi)/4 - (eta*zeta)/4 + (ksi*zeta)/8 - (ksi*zeta^2)/8 + (ksi^2*zeta)/8 - ksi^2/8 - zeta^2/8 - (eta*ksi*zeta)/4 + 1/8, eta/4 + (eta*ksi)/4 - (eta*zeta)/4 - (ksi*zeta)/8 + (ksi*zeta^2)/8 - (ksi^2*zeta)/8 + ksi^2/8 + zeta^2/8 - (eta*ksi*zeta)/4 - 1/8, eta/4 - (eta*ksi)/4 - (eta*zeta)/4 + (ksi*zeta)/8 - (ksi*zeta^2)/8 - (ksi^2*zeta)/8 + ksi^2/8 + zeta^2/8 + (eta*ksi*zeta)/4 - 1/8, eta/4 - (eta*ksi)/4 + (eta*zeta)/4 + (ksi*zeta)/8 + (ksi*zeta^2)/8 - (ksi^2*zeta)/8 - ksi^2/8 - zeta^2/8 - (eta*ksi*zeta)/4 + 1/8, eta/4 + (eta*ksi)/4 + (eta*zeta)/4 - (ksi*zeta)/8 - (ksi*zeta^2)/8 - (ksi^2*zeta)/8 - ksi^2/8 - zeta^2/8 + (eta*ksi*zeta)/4 + 1/8, eta/4 + (eta*ksi)/4 + (eta*zeta)/4 + (ksi*zeta)/8 + (ksi*zeta^2)/8 + (ksi^2*zeta)/8 + ksi^2/8 + zeta^2/8 + (eta*ksi*zeta)/4 - 1/8, eta/4 - (eta*ksi)/4 + (eta*zeta)/4 - (ksi*zeta)/8 - (ksi*zeta^2)/8 + (ksi^2*zeta)/8 + ksi^2/8 + zeta^2/8 - (eta*ksi*zeta)/4 - 1/8,               zeta/4 - (ksi^2*zeta)/4 + ksi^2/4 - 1/4, (eta*zeta)/2 - (eta*ksi)/2 - eta/2 + (eta*ksi*zeta)/2,               (ksi^2*zeta)/4 - zeta/4 - ksi^2/4 + 1/4, (eta*ksi)/2 - eta/2 + (eta*zeta)/2 - (eta*ksi*zeta)/2,                 ksi/4 - (ksi*zeta^2)/4 + zeta^2/4 - 1/4,                 (ksi*zeta^2)/4 - ksi/4 + zeta^2/4 - 1/4,                   ksi/4 - (ksi*zeta^2)/4 - zeta^2/4 + 1/4,                 (ksi*zeta^2)/4 - ksi/4 - zeta^2/4 + 1/4,               (ksi^2*zeta)/4 - zeta/4 + ksi^2/4 - 1/4, - eta/2 - (eta*ksi)/2 - (eta*zeta)/2 - (eta*ksi*zeta)/2,                 zeta/4 - (ksi^2*zeta)/4 - ksi^2/4 + 1/4, (eta*ksi)/2 - eta/2 - (eta*zeta)/2 + (eta*ksi*zeta)/2
            zeta/4 - (eta*ksi)/8 - (eta*zeta)/4 - (ksi*zeta)/4 + (eta*ksi^2)/8 + (eta^2*ksi)/8 - eta^2/8 - ksi^2/8 + (eta*ksi*zeta)/4 + 1/8,   zeta/4 + (eta*ksi)/8 - (eta*zeta)/4 + (ksi*zeta)/4 + (eta*ksi^2)/8 - (eta^2*ksi)/8 - eta^2/8 - ksi^2/8 - (eta*ksi*zeta)/4 + 1/8,   zeta/4 - (eta*ksi)/8 + (eta*zeta)/4 + (ksi*zeta)/4 - (eta*ksi^2)/8 - (eta^2*ksi)/8 - eta^2/8 - ksi^2/8 + (eta*ksi*zeta)/4 + 1/8,   zeta/4 + (eta*ksi)/8 + (eta*zeta)/4 - (ksi*zeta)/4 - (eta*ksi^2)/8 + (eta^2*ksi)/8 - eta^2/8 - ksi^2/8 - (eta*ksi*zeta)/4 + 1/8,   zeta/4 + (eta*ksi)/8 - (eta*zeta)/4 - (ksi*zeta)/4 - (eta*ksi^2)/8 - (eta^2*ksi)/8 + eta^2/8 + ksi^2/8 + (eta*ksi*zeta)/4 - 1/8,   zeta/4 - (eta*ksi)/8 - (eta*zeta)/4 + (ksi*zeta)/4 - (eta*ksi^2)/8 + (eta^2*ksi)/8 + eta^2/8 + ksi^2/8 - (eta*ksi*zeta)/4 - 1/8,   zeta/4 + (eta*ksi)/8 + (eta*zeta)/4 + (ksi*zeta)/4 + (eta*ksi^2)/8 + (eta^2*ksi)/8 + eta^2/8 + ksi^2/8 + (eta*ksi*zeta)/4 - 1/8,   zeta/4 - (eta*ksi)/8 + (eta*zeta)/4 - (ksi*zeta)/4 + (eta*ksi^2)/8 - (eta^2*ksi)/8 + eta^2/8 + ksi^2/8 - (eta*ksi*zeta)/4 - 1/8,                 eta/4 - (eta*ksi^2)/4 + ksi^2/4 - 1/4,                 (eta^2*ksi)/4 - ksi/4 + eta^2/4 - 1/4,                 (eta*ksi^2)/4 - eta/4 + ksi^2/4 - 1/4,                 ksi/4 - (eta^2*ksi)/4 + eta^2/4 - 1/4, (eta*zeta)/2 - zeta/2 + (ksi*zeta)/2 - (eta*ksi*zeta)/2, (eta*zeta)/2 - zeta/2 - (ksi*zeta)/2 + (eta*ksi*zeta)/2, - zeta/2 - (eta*zeta)/2 - (ksi*zeta)/2 - (eta*ksi*zeta)/2, (ksi*zeta)/2 - (eta*zeta)/2 - zeta/2 + (eta*ksi*zeta)/2,                 (eta*ksi^2)/4 - eta/4 - ksi^2/4 + 1/4,                   ksi/4 - (eta^2*ksi)/4 - eta^2/4 + 1/4,                   eta/4 - (eta*ksi^2)/4 - ksi^2/4 + 1/4,                 (eta^2*ksi)/4 - ksi/4 - eta^2/4 + 1/4];
        
        % Derivadas de x,y, respecto de ksi, eta
        jac = dN*nodesEle;
        
        dNxyz = jac\dN;
        
        B = zeros(size(C,2),nDofNod*nNodEle);
        B(1,1:3:(length(B)-2)) = dNxyz(1,:);
        B(2,2:3:(length(B)-1)) = dNxyz(2,:);
        B(3,3:3:length(B))     = dNxyz(3,:);
        B(4,1:3:(length(B)-2)) = dNxyz(2,:);
        B(4,2:3:(length(B)-1)) = dNxyz(1,:);
        B(5,2:3:(length(B)-1)) = dNxyz(3,:);
        B(5,3:3:length(B))     = dNxyz(2,:);
        B(6,1:3:(length(B)-2)) = dNxyz(3,:);
        B(6,3:3:length(B))     = dNxyz(1,:);
        
        eleDofs = nodeDofs(elements(iele,:),:);
        eleDofs = reshape(eleDofs',[],1);
        StressGauss8(iele,ipg,:) = C*B*D(eleDofs)-C*epsilon;
    end
end
%% Extrapolo a Nodos las tensiones a nodos
Nodos=[-1    -1   -1
        1    -1   -1
        1     1   -1
       -1     1   -1
       -1    -1    1
        1    -1    1
        1     1    1
       -1     1    1
        0    -1   -1
        1     0   -1
        0     1   -1
       -1     0   -1
       -1    -1    0
        1    -1    0
        1     1    0
       -1     1    0
        0    -1    1
        1     0    1
        0     1    1
       -1     0    1]*(3)^0.5;
Func= @(ksi,eta,zeta) [ (ksi.*eta)/8 - eta/8 - zeta/8 - ksi/8 + (ksi.*zeta)/8 + (eta.*zeta)/8 - (ksi.*eta.*zeta)/8 + 1/8, ksi/8 - eta/8 - zeta/8 - (ksi.*eta)/8 - (ksi.*zeta)/8 + (eta.*zeta)/8 + (ksi.*eta.*zeta)/8 + 1/8, ksi/8 + eta/8 - zeta/8 + (ksi.*eta)/8 - (ksi.*zeta)/8 - (eta.*zeta)/8 - (ksi.*eta.*zeta)/8 + 1/8, eta/8 - ksi/8 - zeta/8 - (ksi.*eta)/8 + (ksi.*zeta)/8 - (eta.*zeta)/8 + (ksi.*eta.*zeta)/8 + 1/8, zeta/8 - eta/8 - ksi/8 + (ksi.*eta)/8 - (ksi.*zeta)/8 - (eta.*zeta)/8 + (ksi.*eta.*zeta)/8 + 1/8, ksi/8 - eta/8 + zeta/8 - (ksi.*eta)/8 + (ksi.*zeta)/8 - (eta.*zeta)/8 - (ksi.*eta.*zeta)/8 + 1/8, ksi/8 + eta/8 + zeta/8 + (ksi.*eta)/8 + (ksi.*zeta)/8 + (eta.*zeta)/8 + (ksi.*eta.*zeta)/8 + 1/8, eta/8 - ksi/8 + zeta/8 - (ksi.*eta)/8 - (ksi.*zeta)/8 + (eta.*zeta)/8 - (ksi.*eta.*zeta)/8 + 1/8];
N= Func(Nodos(:,1),Nodos(:,2),Nodos(:,3));
stressExt8= zeros(nel,nNodEle,6);
AuxStressGauss8=permute(StressGauss8,[2 3 1]);
for iele = 1:nel
    stressExt8(iele,:,:) = N*AuxStressGauss8(:,:,iele);
end


%% Recuperación de tensiones en los Centros
StressCentro = zeros(nel,nNodEle,6);
for iele = 1:nel
    nodesEle = nodes(elements(iele,:),:);
    if max(nodes(elements(iele,:),3))<= Espesor1
        C=C1;
        epsilon=alfa1*T*[1;1;1;0;0;0];
    else
        C=C2;
        epsilon=alfa2*T*[1;1;1;0;0;0];
    end
    for inode = 1:nNodEle
        % Punto de Gauss
        ksi = 0;
        eta = 0;
        zeta = 0;
        % Derivadas de las funciones de forma respecto de ksi, eta
        dN=[ ksi/4 - (eta*ksi)/4 - (eta*zeta)/8 - (ksi*zeta)/4 + (eta*zeta^2)/8 + (eta^2*zeta)/8 - eta^2/8 - zeta^2/8 + (eta*ksi*zeta)/4 + 1/8, ksi/4 - (eta*ksi)/4 + (eta*zeta)/8 - (ksi*zeta)/4 - (eta*zeta^2)/8 - (eta^2*zeta)/8 + eta^2/8 + zeta^2/8 + (eta*ksi*zeta)/4 - 1/8, ksi/4 + (eta*ksi)/4 - (eta*zeta)/8 - (ksi*zeta)/4 + (eta*zeta^2)/8 - (eta^2*zeta)/8 + eta^2/8 + zeta^2/8 - (eta*ksi*zeta)/4 - 1/8, ksi/4 + (eta*ksi)/4 + (eta*zeta)/8 - (ksi*zeta)/4 - (eta*zeta^2)/8 + (eta^2*zeta)/8 - eta^2/8 - zeta^2/8 - (eta*ksi*zeta)/4 + 1/8, ksi/4 - (eta*ksi)/4 + (eta*zeta)/8 + (ksi*zeta)/4 + (eta*zeta^2)/8 - (eta^2*zeta)/8 - eta^2/8 - zeta^2/8 - (eta*ksi*zeta)/4 + 1/8, ksi/4 - (eta*ksi)/4 - (eta*zeta)/8 + (ksi*zeta)/4 - (eta*zeta^2)/8 + (eta^2*zeta)/8 + eta^2/8 + zeta^2/8 - (eta*ksi*zeta)/4 - 1/8, ksi/4 + (eta*ksi)/4 + (eta*zeta)/8 + (ksi*zeta)/4 + (eta*zeta^2)/8 + (eta^2*zeta)/8 + eta^2/8 + zeta^2/8 + (eta*ksi*zeta)/4 - 1/8, ksi/4 + (eta*ksi)/4 - (eta*zeta)/8 + (ksi*zeta)/4 - (eta*zeta^2)/8 - (eta^2*zeta)/8 - eta^2/8 - zeta^2/8 + (eta*ksi*zeta)/4 + 1/8, (eta*ksi)/2 - ksi/2 + (ksi*zeta)/2 - (eta*ksi*zeta)/2,               (eta^2*zeta)/4 - zeta/4 - eta^2/4 + 1/4, (ksi*zeta)/2 - (eta*ksi)/2 - ksi/2 + (eta*ksi*zeta)/2,               zeta/4 - (eta^2*zeta)/4 + eta^2/4 - 1/4,                 eta/4 - (eta*zeta^2)/4 + zeta^2/4 - 1/4,                 (eta*zeta^2)/4 - eta/4 - zeta^2/4 + 1/4,                   eta/4 - (eta*zeta^2)/4 - zeta^2/4 + 1/4,                 (eta*zeta^2)/4 - eta/4 + zeta^2/4 - 1/4, (eta*ksi)/2 - ksi/2 - (ksi*zeta)/2 + (eta*ksi*zeta)/2,                 zeta/4 - (eta^2*zeta)/4 - eta^2/4 + 1/4, - ksi/2 - (eta*ksi)/2 - (ksi*zeta)/2 - (eta*ksi*zeta)/2,               (eta^2*zeta)/4 - zeta/4 + eta^2/4 - 1/4
            eta/4 - (eta*ksi)/4 - (eta*zeta)/4 - (ksi*zeta)/8 + (ksi*zeta^2)/8 + (ksi^2*zeta)/8 - ksi^2/8 - zeta^2/8 + (eta*ksi*zeta)/4 + 1/8, eta/4 + (eta*ksi)/4 - (eta*zeta)/4 + (ksi*zeta)/8 - (ksi*zeta^2)/8 + (ksi^2*zeta)/8 - ksi^2/8 - zeta^2/8 - (eta*ksi*zeta)/4 + 1/8, eta/4 + (eta*ksi)/4 - (eta*zeta)/4 - (ksi*zeta)/8 + (ksi*zeta^2)/8 - (ksi^2*zeta)/8 + ksi^2/8 + zeta^2/8 - (eta*ksi*zeta)/4 - 1/8, eta/4 - (eta*ksi)/4 - (eta*zeta)/4 + (ksi*zeta)/8 - (ksi*zeta^2)/8 - (ksi^2*zeta)/8 + ksi^2/8 + zeta^2/8 + (eta*ksi*zeta)/4 - 1/8, eta/4 - (eta*ksi)/4 + (eta*zeta)/4 + (ksi*zeta)/8 + (ksi*zeta^2)/8 - (ksi^2*zeta)/8 - ksi^2/8 - zeta^2/8 - (eta*ksi*zeta)/4 + 1/8, eta/4 + (eta*ksi)/4 + (eta*zeta)/4 - (ksi*zeta)/8 - (ksi*zeta^2)/8 - (ksi^2*zeta)/8 - ksi^2/8 - zeta^2/8 + (eta*ksi*zeta)/4 + 1/8, eta/4 + (eta*ksi)/4 + (eta*zeta)/4 + (ksi*zeta)/8 + (ksi*zeta^2)/8 + (ksi^2*zeta)/8 + ksi^2/8 + zeta^2/8 + (eta*ksi*zeta)/4 - 1/8, eta/4 - (eta*ksi)/4 + (eta*zeta)/4 - (ksi*zeta)/8 - (ksi*zeta^2)/8 + (ksi^2*zeta)/8 + ksi^2/8 + zeta^2/8 - (eta*ksi*zeta)/4 - 1/8,               zeta/4 - (ksi^2*zeta)/4 + ksi^2/4 - 1/4, (eta*zeta)/2 - (eta*ksi)/2 - eta/2 + (eta*ksi*zeta)/2,               (ksi^2*zeta)/4 - zeta/4 - ksi^2/4 + 1/4, (eta*ksi)/2 - eta/2 + (eta*zeta)/2 - (eta*ksi*zeta)/2,                 ksi/4 - (ksi*zeta^2)/4 + zeta^2/4 - 1/4,                 (ksi*zeta^2)/4 - ksi/4 + zeta^2/4 - 1/4,                   ksi/4 - (ksi*zeta^2)/4 - zeta^2/4 + 1/4,                 (ksi*zeta^2)/4 - ksi/4 - zeta^2/4 + 1/4,               (ksi^2*zeta)/4 - zeta/4 + ksi^2/4 - 1/4, - eta/2 - (eta*ksi)/2 - (eta*zeta)/2 - (eta*ksi*zeta)/2,                 zeta/4 - (ksi^2*zeta)/4 - ksi^2/4 + 1/4, (eta*ksi)/2 - eta/2 - (eta*zeta)/2 + (eta*ksi*zeta)/2
            zeta/4 - (eta*ksi)/8 - (eta*zeta)/4 - (ksi*zeta)/4 + (eta*ksi^2)/8 + (eta^2*ksi)/8 - eta^2/8 - ksi^2/8 + (eta*ksi*zeta)/4 + 1/8,   zeta/4 + (eta*ksi)/8 - (eta*zeta)/4 + (ksi*zeta)/4 + (eta*ksi^2)/8 - (eta^2*ksi)/8 - eta^2/8 - ksi^2/8 - (eta*ksi*zeta)/4 + 1/8,   zeta/4 - (eta*ksi)/8 + (eta*zeta)/4 + (ksi*zeta)/4 - (eta*ksi^2)/8 - (eta^2*ksi)/8 - eta^2/8 - ksi^2/8 + (eta*ksi*zeta)/4 + 1/8,   zeta/4 + (eta*ksi)/8 + (eta*zeta)/4 - (ksi*zeta)/4 - (eta*ksi^2)/8 + (eta^2*ksi)/8 - eta^2/8 - ksi^2/8 - (eta*ksi*zeta)/4 + 1/8,   zeta/4 + (eta*ksi)/8 - (eta*zeta)/4 - (ksi*zeta)/4 - (eta*ksi^2)/8 - (eta^2*ksi)/8 + eta^2/8 + ksi^2/8 + (eta*ksi*zeta)/4 - 1/8,   zeta/4 - (eta*ksi)/8 - (eta*zeta)/4 + (ksi*zeta)/4 - (eta*ksi^2)/8 + (eta^2*ksi)/8 + eta^2/8 + ksi^2/8 - (eta*ksi*zeta)/4 - 1/8,   zeta/4 + (eta*ksi)/8 + (eta*zeta)/4 + (ksi*zeta)/4 + (eta*ksi^2)/8 + (eta^2*ksi)/8 + eta^2/8 + ksi^2/8 + (eta*ksi*zeta)/4 - 1/8,   zeta/4 - (eta*ksi)/8 + (eta*zeta)/4 - (ksi*zeta)/4 + (eta*ksi^2)/8 - (eta^2*ksi)/8 + eta^2/8 + ksi^2/8 - (eta*ksi*zeta)/4 - 1/8,                 eta/4 - (eta*ksi^2)/4 + ksi^2/4 - 1/4,                 (eta^2*ksi)/4 - ksi/4 + eta^2/4 - 1/4,                 (eta*ksi^2)/4 - eta/4 + ksi^2/4 - 1/4,                 ksi/4 - (eta^2*ksi)/4 + eta^2/4 - 1/4, (eta*zeta)/2 - zeta/2 + (ksi*zeta)/2 - (eta*ksi*zeta)/2, (eta*zeta)/2 - zeta/2 - (ksi*zeta)/2 + (eta*ksi*zeta)/2, - zeta/2 - (eta*zeta)/2 - (ksi*zeta)/2 - (eta*ksi*zeta)/2, (ksi*zeta)/2 - (eta*zeta)/2 - zeta/2 + (eta*ksi*zeta)/2,                 (eta*ksi^2)/4 - eta/4 - ksi^2/4 + 1/4,                   ksi/4 - (eta^2*ksi)/4 - eta^2/4 + 1/4,                   eta/4 - (eta*ksi^2)/4 - ksi^2/4 + 1/4,                 (eta^2*ksi)/4 - ksi/4 - eta^2/4 + 1/4];
        
        % Derivadas de x,y, respecto de ksi, eta
        jac = dN*nodesEle;
        % Derivadas de las funciones de forma respecto de x,y.
        dNxyz = jac\dN;          % dNxy = inv(jac)*dN
        
        B = zeros(size(C,2),nDofNod*nNodEle);
        B(1,1:3:(length(B)-2)) = dNxyz(1,:);
        B(2,2:3:(length(B)-1)) = dNxyz(2,:);
        B(3,3:3:length(B))     = dNxyz(3,:);
        B(4,1:3:(length(B)-2)) = dNxyz(2,:);
        B(4,2:3:(length(B)-1)) = dNxyz(1,:);
        B(5,2:3:(length(B)-1)) = dNxyz(3,:);
        B(5,3:3:length(B))     = dNxyz(2,:);
        B(6,1:3:(length(B)-2)) = dNxyz(3,:);
        B(6,3:3:length(B))     = dNxyz(1,:);
        
        eleDofs = nodeDofs(elements(iele,:),:);
        eleDofs = reshape(eleDofs',[],1);
        StressCentro(iele,inode,:) = C*B*D(eleDofs)-C*epsilon;
    end
end





%% Promedio tensiones en los nodos

avgStress8=zeros(nNod,6);
for inode=1:nNod
    [I,J]=find(elements == inode);
    nShare=length(I);
    for ishare=1:nShare
        avgStress8(inode,:)= avgStress8(inode,:)+(squeeze(stressExt8(I(ishare),J(ishare),:))');%+(squeeze(StressCentro(I(ishare),J(ishare),:))'))/2;
    end
    avgStress8(inode,:)= avgStress8(inode,:)/nShare;
end


avgStressCentro=zeros(nNod,6);
for inode=1:nNod
    [I,J]=find(elements == inode);
    nShare=length(I);
    for ishare=1:nShare
        avgStressCentro(inode,:)= avgStressCentro(inode,:)+(squeeze(stressExt8(I(ishare),J(ishare),:))'+(squeeze(StressCentro(I(ishare),J(ishare),:))'))/2;
    end
    avgStressCentro(inode,:)= avgStressCentro(inode,:)/nShare;
end




%% Esta seccion es para plotear no mas
AvgStress8=zeros(nel,20,6);
for i=1:nel
    AvgStress8(i,:,1)=avgStress8(elements(i,:),1);
    AvgStress8(i,:,2)=avgStress8(elements(i,:),2);
    AvgStress8(i,:,3)=avgStress8(elements(i,:),3);
    AvgStress8(i,:,4)=avgStress8(elements(i,:),4);
    AvgStress8(i,:,5)=avgStress8(elements(i,:),5);
    AvgStress8(i,:,6)=avgStress8(elements(i,:),6);
end


AvgStressCentro=zeros(nel,20,6);
for i=1:nel
    AvgStressCentro(i,:,1)=avgStressCentro(elements(i,:),1);
    AvgStressCentro(i,:,2)=avgStressCentro(elements(i,:),2);
    AvgStressCentro(i,:,3)=avgStressCentro(elements(i,:),3);
    AvgStressCentro(i,:,4)=avgStressCentro(elements(i,:),4);
    AvgStressCentro(i,:,5)=avgStressCentro(elements(i,:),5);
    AvgStressCentro(i,:,6)=avgStressCentro(elements(i,:),6);
end


DisplacementModule=sqrt(diag((Dnodes-nodes)*(Dnodes-nodes)'));
MaxDisplacement=max(DisplacementModule);
vonMisses=(0.5*((stressExt8(:,:,1)-stressExt8(:,:,2)).^2+(stressExt8(:,:,2)-stressExt8(:,:,3)).^2+(stressExt8(:,:,3)-stressExt8(:,:,1)).^2+6*(stressExt8(:,:,4).^2+stressExt8(:,:,5).^2+stressExt8(:,:,6).^2))).^0.5;
vonMisses2=(0.5*((StressGauss8(:,:,1)-StressGauss8(:,:,2)).^2+(StressGauss8(:,:,2)-StressGauss8(:,:,3)).^2+(StressGauss8(:,:,3)-StressGauss8(:,:,1)).^2+6*(StressGauss8(:,:,4).^2+StressGauss8(:,:,5).^2+StressGauss8(:,:,6).^2))).^0.5;
vonMisses3=(stressExt27(:,:,1).^2+stressExt27(:,:,2).^2+stressExt27(:,:,3).^2-(stressExt27(:,:,1).*stressExt27(:,:,2)+stressExt27(:,:,2).*stressExt27(:,:,3)+stressExt27(:,:,3).*stressExt27(:,:,1))+3*(stressExt27(:,:,4).^2+stressExt27(:,:,5).^2+stressExt27(:,:,6).^2)).^0.5;
vonMisses4=(StressGauss27(:,:,1).^2+StressGauss27(:,:,2).^2+StressGauss27(:,:,3).^2-(StressGauss27(:,:,1).*StressGauss27(:,:,2)+StressGauss27(:,:,2).*StressGauss27(:,:,3)+StressGauss27(:,:,3).*StressGauss27(:,:,1))+3*(StressGauss27(:,:,4).^2+StressGauss27(:,:,5).^2+StressGauss27(:,:,6).^2)).^0.5;

VonMisses=zeros(nNod,6);
for inode=1:nNod
    [I,J]=find(elements == inode);
    nShare=length(I);
    for ishare=1:nShare
        VonMisses(inode,:)= VonMisses(inode,:)+(squeeze(vonMisses(I(ishare),J(ishare),:))');%+(squeeze(StressCentro(I(ishare),J(ishare),:))'))/2;
    end
    VonMisses(inode,:)= VonMisses(inode,:)/nShare;
end
AvgVonMisses=zeros(nel,20);
for i=1:nel
    AvgVonMisses(i,:)=VonMisses(elements(i,:));
end

if Espesor1<2
elemUp=(1:(ElementosLargo*ElementosH))+(ElementosLargo*ElementosH)*Espesor1/Espesor*ElementosEspesor;
elemDown=(1:(ElementosLargo*ElementosH))+(ElementosLargo*ElementosH)*(Espesor1/Espesor*ElementosEspesor-1);
nodesDown=[5 6 7 8 17 18 19 20];
nodesUp=[1 2 3 4 9 10 11 12];
avgStress8=AvgStress8;
avgStress8(elemUp,nodesUp,:)=stressExt8(elemUp,nodesUp,:);
avgStress8(elemDown,nodesDown,:)=stressExt8(elemDown,nodesDown,:);
AvgVonMisses(elemUp,nodesUp,:)=vonMisses(elemUp,nodesUp,:);
AvgVonMisses(elemDown,nodesDown,:)=vonMisses(elemDown,nodesDown,:);
VonMisses(elemUp,nodesUp,:)=vonMisses(elemUp,nodesUp,:);
VonMisses(elemDown,nodesDown,:)=vonMisses(elemDown,nodesDown,:);
end


xyzgauss=[-1    -1   -1
        1    -1   -1
        1     1   -1
       -1     1   -1
       -1    -1    1
        1    -1    1
        1     1    1
       -1     1    1];
upg=[xyzgauss]/3^0.5; % Ubicaciones puntos de Gauss
wpg=ones(8,1);
npg = size(upg,1); % Numero de puntos de Gauss
invC1=C1\eye(6);
invC2=C2\eye(6);
Eta_elem=zeros(nel,1);
E2_elem=zeros(nel,1);
U2_elem=zeros(nel,1);
for iele = 1:nel
    nodesEle = nodes(elements(iele,:),:);
    if max(nodes(elements(iele,:),3))<= Espesor1
        invC=invC1;
    else
        invC=invC2;
    end
    for ipg = 1:npg
        % Punto de Gauss
        ksi = upg(ipg,1);
        eta = upg(ipg,2);
        zeta = upg(ipg,3);
        % Derivadas de las funciones de forma respecto de ksi, eta
        %         [N dN]=FuncionFormaH20(ksi,eta,zeta);
      dN=[ ksi/4 - (eta*ksi)/4 - (eta*zeta)/8 - (ksi*zeta)/4 + (eta*zeta^2)/8 + (eta^2*zeta)/8 - eta^2/8 - zeta^2/8 + (eta*ksi*zeta)/4 + 1/8, ksi/4 - (eta*ksi)/4 + (eta*zeta)/8 - (ksi*zeta)/4 - (eta*zeta^2)/8 - (eta^2*zeta)/8 + eta^2/8 + zeta^2/8 + (eta*ksi*zeta)/4 - 1/8, ksi/4 + (eta*ksi)/4 - (eta*zeta)/8 - (ksi*zeta)/4 + (eta*zeta^2)/8 - (eta^2*zeta)/8 + eta^2/8 + zeta^2/8 - (eta*ksi*zeta)/4 - 1/8, ksi/4 + (eta*ksi)/4 + (eta*zeta)/8 - (ksi*zeta)/4 - (eta*zeta^2)/8 + (eta^2*zeta)/8 - eta^2/8 - zeta^2/8 - (eta*ksi*zeta)/4 + 1/8, ksi/4 - (eta*ksi)/4 + (eta*zeta)/8 + (ksi*zeta)/4 + (eta*zeta^2)/8 - (eta^2*zeta)/8 - eta^2/8 - zeta^2/8 - (eta*ksi*zeta)/4 + 1/8, ksi/4 - (eta*ksi)/4 - (eta*zeta)/8 + (ksi*zeta)/4 - (eta*zeta^2)/8 + (eta^2*zeta)/8 + eta^2/8 + zeta^2/8 - (eta*ksi*zeta)/4 - 1/8, ksi/4 + (eta*ksi)/4 + (eta*zeta)/8 + (ksi*zeta)/4 + (eta*zeta^2)/8 + (eta^2*zeta)/8 + eta^2/8 + zeta^2/8 + (eta*ksi*zeta)/4 - 1/8, ksi/4 + (eta*ksi)/4 - (eta*zeta)/8 + (ksi*zeta)/4 - (eta*zeta^2)/8 - (eta^2*zeta)/8 - eta^2/8 - zeta^2/8 + (eta*ksi*zeta)/4 + 1/8, (eta*ksi)/2 - ksi/2 + (ksi*zeta)/2 - (eta*ksi*zeta)/2,               (eta^2*zeta)/4 - zeta/4 - eta^2/4 + 1/4, (ksi*zeta)/2 - (eta*ksi)/2 - ksi/2 + (eta*ksi*zeta)/2,               zeta/4 - (eta^2*zeta)/4 + eta^2/4 - 1/4,                 eta/4 - (eta*zeta^2)/4 + zeta^2/4 - 1/4,                 (eta*zeta^2)/4 - eta/4 - zeta^2/4 + 1/4,                   eta/4 - (eta*zeta^2)/4 - zeta^2/4 + 1/4,                 (eta*zeta^2)/4 - eta/4 + zeta^2/4 - 1/4, (eta*ksi)/2 - ksi/2 - (ksi*zeta)/2 + (eta*ksi*zeta)/2,                 zeta/4 - (eta^2*zeta)/4 - eta^2/4 + 1/4, - ksi/2 - (eta*ksi)/2 - (ksi*zeta)/2 - (eta*ksi*zeta)/2,               (eta^2*zeta)/4 - zeta/4 + eta^2/4 - 1/4
            eta/4 - (eta*ksi)/4 - (eta*zeta)/4 - (ksi*zeta)/8 + (ksi*zeta^2)/8 + (ksi^2*zeta)/8 - ksi^2/8 - zeta^2/8 + (eta*ksi*zeta)/4 + 1/8, eta/4 + (eta*ksi)/4 - (eta*zeta)/4 + (ksi*zeta)/8 - (ksi*zeta^2)/8 + (ksi^2*zeta)/8 - ksi^2/8 - zeta^2/8 - (eta*ksi*zeta)/4 + 1/8, eta/4 + (eta*ksi)/4 - (eta*zeta)/4 - (ksi*zeta)/8 + (ksi*zeta^2)/8 - (ksi^2*zeta)/8 + ksi^2/8 + zeta^2/8 - (eta*ksi*zeta)/4 - 1/8, eta/4 - (eta*ksi)/4 - (eta*zeta)/4 + (ksi*zeta)/8 - (ksi*zeta^2)/8 - (ksi^2*zeta)/8 + ksi^2/8 + zeta^2/8 + (eta*ksi*zeta)/4 - 1/8, eta/4 - (eta*ksi)/4 + (eta*zeta)/4 + (ksi*zeta)/8 + (ksi*zeta^2)/8 - (ksi^2*zeta)/8 - ksi^2/8 - zeta^2/8 - (eta*ksi*zeta)/4 + 1/8, eta/4 + (eta*ksi)/4 + (eta*zeta)/4 - (ksi*zeta)/8 - (ksi*zeta^2)/8 - (ksi^2*zeta)/8 - ksi^2/8 - zeta^2/8 + (eta*ksi*zeta)/4 + 1/8, eta/4 + (eta*ksi)/4 + (eta*zeta)/4 + (ksi*zeta)/8 + (ksi*zeta^2)/8 + (ksi^2*zeta)/8 + ksi^2/8 + zeta^2/8 + (eta*ksi*zeta)/4 - 1/8, eta/4 - (eta*ksi)/4 + (eta*zeta)/4 - (ksi*zeta)/8 - (ksi*zeta^2)/8 + (ksi^2*zeta)/8 + ksi^2/8 + zeta^2/8 - (eta*ksi*zeta)/4 - 1/8,               zeta/4 - (ksi^2*zeta)/4 + ksi^2/4 - 1/4, (eta*zeta)/2 - (eta*ksi)/2 - eta/2 + (eta*ksi*zeta)/2,               (ksi^2*zeta)/4 - zeta/4 - ksi^2/4 + 1/4, (eta*ksi)/2 - eta/2 + (eta*zeta)/2 - (eta*ksi*zeta)/2,                 ksi/4 - (ksi*zeta^2)/4 + zeta^2/4 - 1/4,                 (ksi*zeta^2)/4 - ksi/4 + zeta^2/4 - 1/4,                   ksi/4 - (ksi*zeta^2)/4 - zeta^2/4 + 1/4,                 (ksi*zeta^2)/4 - ksi/4 - zeta^2/4 + 1/4,               (ksi^2*zeta)/4 - zeta/4 + ksi^2/4 - 1/4, - eta/2 - (eta*ksi)/2 - (eta*zeta)/2 - (eta*ksi*zeta)/2,                 zeta/4 - (ksi^2*zeta)/4 - ksi^2/4 + 1/4, (eta*ksi)/2 - eta/2 - (eta*zeta)/2 + (eta*ksi*zeta)/2
            zeta/4 - (eta*ksi)/8 - (eta*zeta)/4 - (ksi*zeta)/4 + (eta*ksi^2)/8 + (eta^2*ksi)/8 - eta^2/8 - ksi^2/8 + (eta*ksi*zeta)/4 + 1/8,   zeta/4 + (eta*ksi)/8 - (eta*zeta)/4 + (ksi*zeta)/4 + (eta*ksi^2)/8 - (eta^2*ksi)/8 - eta^2/8 - ksi^2/8 - (eta*ksi*zeta)/4 + 1/8,   zeta/4 - (eta*ksi)/8 + (eta*zeta)/4 + (ksi*zeta)/4 - (eta*ksi^2)/8 - (eta^2*ksi)/8 - eta^2/8 - ksi^2/8 + (eta*ksi*zeta)/4 + 1/8,   zeta/4 + (eta*ksi)/8 + (eta*zeta)/4 - (ksi*zeta)/4 - (eta*ksi^2)/8 + (eta^2*ksi)/8 - eta^2/8 - ksi^2/8 - (eta*ksi*zeta)/4 + 1/8,   zeta/4 + (eta*ksi)/8 - (eta*zeta)/4 - (ksi*zeta)/4 - (eta*ksi^2)/8 - (eta^2*ksi)/8 + eta^2/8 + ksi^2/8 + (eta*ksi*zeta)/4 - 1/8,   zeta/4 - (eta*ksi)/8 - (eta*zeta)/4 + (ksi*zeta)/4 - (eta*ksi^2)/8 + (eta^2*ksi)/8 + eta^2/8 + ksi^2/8 - (eta*ksi*zeta)/4 - 1/8,   zeta/4 + (eta*ksi)/8 + (eta*zeta)/4 + (ksi*zeta)/4 + (eta*ksi^2)/8 + (eta^2*ksi)/8 + eta^2/8 + ksi^2/8 + (eta*ksi*zeta)/4 - 1/8,   zeta/4 - (eta*ksi)/8 + (eta*zeta)/4 - (ksi*zeta)/4 + (eta*ksi^2)/8 - (eta^2*ksi)/8 + eta^2/8 + ksi^2/8 - (eta*ksi*zeta)/4 - 1/8,                 eta/4 - (eta*ksi^2)/4 + ksi^2/4 - 1/4,                 (eta^2*ksi)/4 - ksi/4 + eta^2/4 - 1/4,                 (eta*ksi^2)/4 - eta/4 + ksi^2/4 - 1/4,                 ksi/4 - (eta^2*ksi)/4 + eta^2/4 - 1/4, (eta*zeta)/2 - zeta/2 + (ksi*zeta)/2 - (eta*ksi*zeta)/2, (eta*zeta)/2 - zeta/2 - (ksi*zeta)/2 + (eta*ksi*zeta)/2, - zeta/2 - (eta*zeta)/2 - (ksi*zeta)/2 - (eta*ksi*zeta)/2, (ksi*zeta)/2 - (eta*zeta)/2 - zeta/2 + (eta*ksi*zeta)/2,                 (eta*ksi^2)/4 - eta/4 - ksi^2/4 + 1/4,                   ksi/4 - (eta^2*ksi)/4 - eta^2/4 + 1/4,                   eta/4 - (eta*ksi^2)/4 - ksi^2/4 + 1/4,                 (eta^2*ksi)/4 - ksi/4 - eta^2/4 + 1/4];
        % Derivadas de x,y, respecto de ksi, eta
        jac = dN*nodesEle;
        
        dNxyz = jac\dN;
        
        N=[ (ksi.^2.*eta.*zeta)/8 - (ksi.^2.*eta)/8 - (ksi.^2.*zeta)/8 + ksi.^2/8 + (ksi.*eta.^2.*zeta)/8 - (ksi.*eta.^2)/8 + (ksi.*eta.*zeta.^2)/8 - (ksi.*eta.*zeta)/8 - (ksi.*zeta.^2)/8 + ksi/8 - (eta.^2.*zeta)/8 + eta.^2/8 - (eta.*zeta.^2)/8 + eta/8 + zeta.^2/8 + zeta/8 - 1/4, (ksi.^2.*eta.*zeta)/8 - (ksi.^2.*eta)/8 - (ksi.^2.*zeta)/8 + ksi.^2/8 - (ksi.*eta.^2.*zeta)/8 + (ksi.*eta.^2)/8 - (ksi.*eta.*zeta.^2)/8 + (ksi.*eta.*zeta)/8 + (ksi.*zeta.^2)/8 - ksi/8 - (eta.^2.*zeta)/8 + eta.^2/8 - (eta.*zeta.^2)/8 + eta/8 + zeta.^2/8 + zeta/8 - 1/4, - (ksi.^2.*eta.*zeta)/8 + (ksi.^2.*eta)/8 - (ksi.^2.*zeta)/8 + ksi.^2/8 - (ksi.*eta.^2.*zeta)/8 + (ksi.*eta.^2)/8 + (ksi.*eta.*zeta.^2)/8 - (ksi.*eta.*zeta)/8 + (ksi.*zeta.^2)/8 - ksi/8 - (eta.^2.*zeta)/8 + eta.^2/8 + (eta.*zeta.^2)/8 - eta/8 + zeta.^2/8 + zeta/8 - 1/4, - (ksi.^2.*eta.*zeta)/8 + (ksi.^2.*eta)/8 - (ksi.^2.*zeta)/8 + ksi.^2/8 + (ksi.*eta.^2.*zeta)/8 - (ksi.*eta.^2)/8 - (ksi.*eta.*zeta.^2)/8 + (ksi.*eta.*zeta)/8 - (ksi.*zeta.^2)/8 + ksi/8 - (eta.^2.*zeta)/8 + eta.^2/8 + (eta.*zeta.^2)/8 - eta/8 + zeta.^2/8 + zeta/8 - 1/4, - (ksi.^2.*eta.*zeta)/8 - (ksi.^2.*eta)/8 + (ksi.^2.*zeta)/8 + ksi.^2/8 - (ksi.*eta.^2.*zeta)/8 - (ksi.*eta.^2)/8 + (ksi.*eta.*zeta.^2)/8 + (ksi.*eta.*zeta)/8 - (ksi.*zeta.^2)/8 + ksi/8 + (eta.^2.*zeta)/8 + eta.^2/8 - (eta.*zeta.^2)/8 + eta/8 + zeta.^2/8 - zeta/8 - 1/4, - (ksi.^2.*eta.*zeta)/8 - (ksi.^2.*eta)/8 + (ksi.^2.*zeta)/8 + ksi.^2/8 + (ksi.*eta.^2.*zeta)/8 + (ksi.*eta.^2)/8 - (ksi.*eta.*zeta.^2)/8 - (ksi.*eta.*zeta)/8 + (ksi.*zeta.^2)/8 - ksi/8 + (eta.^2.*zeta)/8 + eta.^2/8 - (eta.*zeta.^2)/8 + eta/8 + zeta.^2/8 - zeta/8 - 1/4, (ksi.^2.*eta.*zeta)/8 + (ksi.^2.*eta)/8 + (ksi.^2.*zeta)/8 + ksi.^2/8 + (ksi.*eta.^2.*zeta)/8 + (ksi.*eta.^2)/8 + (ksi.*eta.*zeta.^2)/8 + (ksi.*eta.*zeta)/8 + (ksi.*zeta.^2)/8 - ksi/8 + (eta.^2.*zeta)/8 + eta.^2/8 + (eta.*zeta.^2)/8 - eta/8 + zeta.^2/8 - zeta/8 - 1/4, (ksi.^2.*eta.*zeta)/8 + (ksi.^2.*eta)/8 + (ksi.^2.*zeta)/8 + ksi.^2/8 - (ksi.*eta.^2.*zeta)/8 - (ksi.*eta.^2)/8 - (ksi.*eta.*zeta.^2)/8 - (ksi.*eta.*zeta)/8 - (ksi.*zeta.^2)/8 + ksi/8 + (eta.^2.*zeta)/8 + eta.^2/8 + (eta.*zeta.^2)/8 - eta/8 + zeta.^2/8 - zeta/8 - 1/4, (eta.*zeta)/4 - zeta/4 - eta/4 + (ksi.^2.*eta)/4 + (ksi.^2.*zeta)/4 - ksi.^2/4 - (ksi.^2.*eta.*zeta)/4 + 1/4, ksi/4 - zeta/4 - (ksi.*zeta)/4 - (ksi.*eta.^2)/4 + (eta.^2.*zeta)/4 - eta.^2/4 + (ksi.*eta.^2.*zeta)/4 + 1/4, eta/4 - zeta/4 - (eta.*zeta)/4 - (ksi.^2.*eta)/4 + (ksi.^2.*zeta)/4 - ksi.^2/4 + (ksi.^2.*eta.*zeta)/4 + 1/4, (ksi.*zeta)/4 - zeta/4 - ksi/4 + (ksi.*eta.^2)/4 + (eta.^2.*zeta)/4 - eta.^2/4 - (ksi.*eta.^2.*zeta)/4 + 1/4, (ksi.*eta)/4 - eta/4 - ksi/4 + (ksi.*zeta.^2)/4 + (eta.*zeta.^2)/4 - zeta.^2/4 - (ksi.*eta.*zeta.^2)/4 + 1/4, ksi/4 - eta/4 - (ksi.*eta)/4 - (ksi.*zeta.^2)/4 + (eta.*zeta.^2)/4 - zeta.^2/4 + (ksi.*eta.*zeta.^2)/4 + 1/4, ksi/4 + eta/4 + (ksi.*eta)/4 - (ksi.*zeta.^2)/4 - (eta.*zeta.^2)/4 - zeta.^2/4 - (ksi.*eta.*zeta.^2)/4 + 1/4, eta/4 - ksi/4 - (ksi.*eta)/4 + (ksi.*zeta.^2)/4 - (eta.*zeta.^2)/4 - zeta.^2/4 + (ksi.*eta.*zeta.^2)/4 + 1/4, zeta/4 - eta/4 - (eta.*zeta)/4 + (ksi.^2.*eta)/4 - (ksi.^2.*zeta)/4 - ksi.^2/4 + (ksi.^2.*eta.*zeta)/4 + 1/4, ksi/4 + zeta/4 + (ksi.*zeta)/4 - (ksi.*eta.^2)/4 - (eta.^2.*zeta)/4 - eta.^2/4 - (ksi.*eta.^2.*zeta)/4 + 1/4, eta/4 + zeta/4 + (eta.*zeta)/4 - (ksi.^2.*eta)/4 - (ksi.^2.*zeta)/4 - ksi.^2/4 - (ksi.^2.*eta.*zeta)/4 + 1/4, zeta/4 - ksi/4 - (ksi.*zeta)/4 + (ksi.*eta.^2)/4 - (eta.^2.*zeta)/4 - eta.^2/4 + (ksi.*eta.^2.*zeta)/4 + 1/4];
        ElemStress=squeeze(StressGauss8(iele,ipg,:));
        SmothedStress=(N*squeeze(avgStress8(iele,:,:)))';
        
        E2_elem(iele)=E2_elem(iele)+(SmothedStress - ElemStress)'*invC*(SmothedStress - ElemStress)*wpg(ipg)*det(jac);
        
        U2_elem(iele)=U2_elem(iele)+ElemStress'*invC*ElemStress*wpg(ipg)*det(jac);
        
    end
        
       Eta_elem(iele)=sqrt( E2_elem(iele)/(E2_elem(iele)+U2_elem(iele)));
end
    
Eta_global=sqrt( sum(E2_elem(iele))/(sum(E2_elem(iele))+sum(U2_elem(iele))));


%
%
% Analitica 1 material
%%Deflexion Euler-Bernoulli
% x=linspace(0,Largo,ElementosLargo+1);
% w=P*x.^2.*(3*Largo-x)/(6*E1*Espesor^3*LadoH/12);
% %%Deflexion Timoshenko
% z=0;
% wt=P/(6*E1*Espesor^3*LadoH/12)*(3*NU1*z^2*(Largo-x)+(4+5*NU1)*(Espesor^2).*x/4+(3*Largo-x).*x.^2);
% z=linspace(-Espesor/2,Espesor/2,ElementosEspesor+1);
% wtend=P/(6*E1*Espesor^3*LadoH/12)*(3*NU1.*z.^2*(Largo-Largo/2)+(4+5*NU1)*(Espesor^2).*Largo/2/4+(3*Largo-Largo/2).*(Largo/2).^2);
% a=find(nodes(:,3)==Espesor/2);
% b=find(nodes(:,2)==LadoH/2);
% c=intersect(a,b);
% figure
% hold on
% plot(nodes(c,1),Dnodes(c,3)-nodes(c,3))
% plot(x,w)
% plot(x,wt)
% %% Chequeo Corte
% % for i=1:11
% % a=find(nodes(:,3)==(1+0.1*(i-1)));
% % %  a=find(nodes(:,1)==0.5*Largo);
% % b=find(nodes(:,2)==4);
% % c=intersect(a,b);
% % % 
% % %     figure 
% % hold on
% % plot(nodes(c,1),avgStress8(c,6))
% % end
% % legend('1','2','3','4','5','6','7','8','9','10')
% 
% 
% 
% a=find(nodes(:,1)==Largo/2);
% b=find(nodes(:,2)==LadoH/2);
% c=intersect(a,b);
% figure
% plot(z,wtend)
% hold on
% plot(nodes(c,3)-Espesor/2,Dnodes(c,3)-nodes(c,3))
% 
% figure
% plot(Dnodes(c,3)-nodes(c,3),nodes(c,3)-Espesor/2)
%  hold on
% plot(wtend,z)
