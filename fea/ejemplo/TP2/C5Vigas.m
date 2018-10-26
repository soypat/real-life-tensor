%% C.24 Vigas Pinned Movil
clear ; close all ; clc
%% Datos del problema
b = 100;       %mm Dimesión del lado de un perfil cuadrado
A = b^2;     %mm^2 Area del perfil
Iz = b^4/12; %mm^4 Momento de inercia con respecto al eje z del prefil
H = -[1000 1500 1700];     %mm Altura de los miembro verticales
L = 1600;     %mm Largo de los miembros horizontales
q = -10;       %N/mm Carga aplicada
E = 210000;  %Mpa supuesto
%% Nodos y elementos
CoordenadasNodos = [0 0
                    L H(1)
                    L 0
                    2*L H(2)
                    2*L 0
                    3*L H(3)
                    3*L 0
                    4*L H(2)
                    4*L 0
                    5*L H(1)
                    5*L 0
                    6*L 0];
NumeroNodos = size(CoordenadasNodos,1);
Elementos = [1 2
             1 3
             2 3
             2 4
             3 4
             3 5
             4 5
             4 6
             5 6
             5 7
             6 7
             6 8
             6 9
             7 9
             8 9
             8 10
             8 11
             9 11
             10 11
             10 12
             11 12];
NumeroElementos = size(Elementos,1);
%% DOF
DofPorNodo = 3; %Posibilidad de moverse en eje x e y únicamente
DofTotales = DofPorNodo*NumeroNodos;
Dof = reshape((1:1:DofTotales)',DofPorNodo,NumeroNodos)';
%% Matriz global
Kglobal = zeros(DofTotales);
Memlong = zeros(NumeroElementos,1);
MemT = zeros(6,6,NumeroElementos);
for e= 1:NumeroElementos
    vec = CoordenadasNodos(Elementos(e,2),:) - CoordenadasNodos(Elementos(e,1),:);
    long = norm(vec);
    Memlong(e) = long;
    vecdirec = vec/long;
    Lambda = [vecdirec 0 ;
         -vecdirec(2) vecdirec(1) 0
         0 0 1];
    T = blkdiag(Lambda,Lambda);
    MemT(:,:,e) = T;
    X = A*E/long;
    Y4 = 2*E*Iz/long;
    Y3 = Y4*2;
    Y2 = Y4*3/long;
    Y1 = Y2*2/long;
    Kbarra = X*[1 -1; -1 1];
    Kviga = [Y1 Y2 -Y1 Y2
             Y2 Y3 -Y2 Y4
            -Y1 -Y2 Y1 -Y2
             Y2 Y4 -Y2 Y3];
    Klocal([1 4],[1 4]) = Kbarra;
    Klocal([2 3 5 6],[2 3 5 6]) = Kviga;
    Kelemento = T'*Klocal*T;
    Dofelemento = [Dof(Elementos(e,1),:) Dof(Elementos(e,2),:)];
    Kglobal(Dofelemento,Dofelemento) = Kglobal(Dofelemento,Dofelemento)+Kelemento;
end
%% BC
Fijos = 1&[1 1 0 zeros(1,30) 0 1 0];
Libres = ~Fijos;
%% Cargas
P = zeros(DofPorNodo,NumeroNodos);
for i= 1:NumeroNodos/2
    impar = i*2-1;
    P(2,impar) = q*L;
end
P(3,[1 12])=q*L^2/12*[1 -1];
P = reshape(P,[],1);
%% Reductor
Pred = P(Libres);
Kglobalred = Kglobal(Libres,Libres);
%% Solver
Dred = Kglobalred\Pred;
D = [0 ; 0 ; Dred(1:end-1) ; 0 ; Dred(end)];
figure; hold on
Draw_Barra(Elementos,CoordenadasNodos,'b')
Dp = reshape(D,DofPorNodo,NumeroNodos)';
Dp = Dp(:,1:2);
CoordenadasFinales = CoordenadasNodos+Dp;
Draw_Barra(Elementos,CoordenadasFinales,'r')
%% Funciones de forma
DPorNodo = reshape(D,DofPorNodo,NumeroNodos)';
SigmaAxial = zeros(NumeroElementos,1);
SigmaFlexSup = zeros(NumeroElementos,30);
SigmaTotalSup = SigmaFlexSup;
SigmaTotalInf = SigmaFlexSup;
Memsub = SigmaFlexSup;
for e = 1:NumeroElementos
    Delemento = [DPorNodo(Elementos(e,1),:) DPorNodo(Elementos(e,2),:)]';
    Dlocal = MemT(:,:,e)*Delemento;
    Ba = [-1 1]/Memlong(e);
    SigmaAxial(e) = E*Ba*Dlocal([1 4]);   
    sub = 0:Memlong(e)/29:Memlong(e);
    Memsub(e,:) = sub;
    N1=@(x) -6/Memlong(e)^2+12*x/Memlong(e)^3;
    N2=@(x) -4/Memlong(e)+6*x/Memlong(e)^2;
    N3=@(x) 6/Memlong(e)^2-12*x/Memlong(e)^3;
    N4=@(x) -2/Memlong(e)+6*x/Memlong(e)^2;
    Bf = [N1(sub); N2(sub); N3(sub); N4(sub)]';
    SigmaFlexSup(e,:) = (b/2*E*Bf*Dlocal([2 3 5 6]))';
    SigmaFlexInf = -SigmaFlexSup;
    SigmaTotalSup(e,:) = SigmaFlexSup(e,:)+SigmaAxial(e);
    SigmaTotalInf(e,:) = SigmaFlexInf(e,:)+SigmaAxial(e);
end
SigmaCompuesto = [(SigmaTotalSup); (SigmaTotalInf)];
SigmaMax = max(max(abs(SigmaCompuesto)));
[i,j] = find(abs(SigmaCompuesto)==SigmaMax);
disp(strcat(['La tensión máxima es de: ',num2str(SigmaCompuesto(i,j)),' MPa']))
if i<22
    disp(strcat(['Se encuentra en el elemento N: ',num2str(i)]))
    disp(strcat(['En la posición: ',num2str(Memsub(j)),' en la parte superior']))
else
    disp(strcat(['Se encuentra en el elemento N: ',num2str(i-21)]))
    disp(strcat(['En la posición: ',num2str(Memsub(j)),' en la parte inferior']))
end