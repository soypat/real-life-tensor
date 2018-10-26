%% C.24 Barras Pin Movil
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
DofPorNodo = 2; %Posibilidad de moverse en eje x e y únicamente
DofTotales = DofPorNodo*NumeroNodos;
Dof = reshape((1:1:DofTotales)',DofPorNodo,NumeroNodos)';
%% Matriz global
Kglobal = zeros(DofTotales);
Memlong = zeros(NumeroElementos,1);
MemT = zeros(2,4,NumeroElementos);
for e = 1:NumeroElementos
    vec = CoordenadasNodos(Elementos(e,2),:) - CoordenadasNodos(Elementos(e,1),:);
    long = norm(vec);
    vecdirec = vec/long;
    Memlong(e) = long;
    T = [vecdirec 0 0
         0 0 vecdirec];
    MemT(:,:,e) = T;
    Klocal = E*A/long*[1 -1; -1 1];
    Kelemento = T'*Klocal*T;
    Dofelemento = [Dof(Elementos(e,1),:) Dof(Elementos(e,2),:)];
    Kglobal(Dofelemento,Dofelemento) = Kglobal(Dofelemento,Dofelemento)+Kelemento;
end
%% BC
Fijos = 1&[1 1 zeros(1,21) 1];
Libres = ~Fijos;
%% Cargas
P = zeros(DofPorNodo,NumeroNodos);
for i= 2:NumeroNodos/2
    impar = i*2-1;
    P(2,impar) = q*L;
end
P = reshape(P,[],1);
%% Reductor
Pred = P(Libres);
Kglobalred = Kglobal(Libres,Libres);
%% Solver
Dred = Kglobalred\Pred;
D = [0 ; 0 ; Dred ; 0];
figure; hold on
Draw_Barra(Elementos,CoordenadasNodos,'b')
CoordenadasFinales = CoordenadasNodos+reshape(D,DofPorNodo,NumeroNodos)';
Draw_Barra(Elementos,CoordenadasFinales,'r')
P = Kglobal*D;
%% Tensiones 1
DPorNodo = reshape(D,DofPorNodo,NumeroNodos)';
Sigma = zeros(NumeroElementos,1);
for e = 1:NumeroElementos
    Delemento = [DPorNodo(Elementos(e,1),:) DPorNodo(Elementos(e,2),:)]';
    Dlocal = MemT(:,:,e)*Delemento;
    B = [-1 1]/Memlong(e);
    Sigma(e) = E*B*Dlocal;
end
SigmaMaxAbs = max(abs(Sigma));
Elemento = find(Sigma==SigmaMaxAbs);
SigmaMax = Sigma(Elemento);