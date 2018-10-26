clear
close all
clc


%% Datos

E = 70e3;     %MPa
NU = 0.3; 
a = 1;       
P = 10e3;      %N

% Nodos
Nodos=[0 0
       1 0
       2 0
       3 0
       4 0
       5 0
       6 0
       7 0
       0 1
       1 1
       2 1
       3 1
       4 1
       5 1
       6 1
       7 1
       0 2
       1 2
       2 2
       3 2
       4 2
       5 2
       6 2
       7 2
       0 3
       1 3
       2 3
       3 3
       4 3
       5 3
       6 3
       7 3
       0 4
       1 4
       2 4
       3 4
       4 4
       5 4
       0 5
       1 5
       2 5
       3 5
       4 5
       5 5]*a;
NroNodos = size(Nodos,1);    

% Elementos
Elementos = [1 2 10 9
             2 3 11 10
             3 4 12 11
             4 5 13 12
             5 6 14 13
             6 7 15 14
             7 8 16 15
             9 10 18 17
             10 11 19 18
             11 12 20 19
             12 13 21 20
             13 14 22 21
             14 15 23 22
             15 16 24 23
             17 18 26 25
             18 19 27 26
             19 20 28 27
             20 21 29 28
             21 22 30 29
             22 23 31 30
             23 24 32 31
             25 26 34 33
             26 27 35 34
             27 28 36 35
             28 29 37 36
             29 30 38 37
             33 34 40 39
             34 35 41 40
             35 36 42 41
             36 37 43 42
             37 38 44 43];
NroEle = size(Elementos,1);
NroNodosEle = 4;  

% Grados de Libertad
GdL = 2;
GdLTot = GdL*NroNodos;
bc = zeros(NroNodos,GdL);
bc(1:6,1:2) = 1;
bc(9,1:2) = 1;
bc(17,1:2) = 1;
bc(25,1:2) = 1;
bc(33,1:2) = 1;
bc(39,1:2) = 1;
bc = logical(reshape(bc',[],1));

% Cargas
R = zeros(NroNodos,GdL);
R(8,2) = -P;
R = reshape(R',[],1);


%% Matriz Constitutiva (plane stress)

C = E/(1 - NU^2)*[ 1.0     NU         0.0
                    NU    1.0         0.0
                   0.0    0.0     (1 - NU)/2 ];

               
%% Gauss
a   = 1/sqrt(3);
% Ubicaciones puntos de Gauss
upg = [ -a  -a
         a  -a
         a   a
        -a   a ];    
% Número de puntos de Gauss
npg = size(upg,1);
wpg = ones(npg,1);


%% Matriz de rigidez
K = zeros(GdLTot);
GdLNodos = reshape(1:GdLTot,GdL,NroNodos)';

for iele = 1:NroEle
    Ke = zeros(GdL*NroNodosEle);
    NodosEle = Nodos(Elementos(iele,:),:);
    
    for ipg = 1:npg
        % Punto de Gauss
        ksi = upg(ipg,1);
        eta = upg(ipg,2);  
        % Derivadas de las funciones de forma respecto de ksi, eta
        dN = 1/4*[-(1-eta)   1-eta    1+eta  -(1+eta)
                  -(1-ksi) -(1+ksi)   1+ksi    1-ksi ];  
        % Derivadas de x,y, respecto de ksi, eta
        jac = dN*NodosEle;                      
        % Derivadas de las funciones de forma respecto de x,y.
        dNxy = jac\dN;          % dNxy = inv(jac)*dN
        
        B = zeros(size(C,2),GdL*NroNodosEle);
        B(1,1:2:7) = dNxy(1,:);
        B(2,2:2:8) = dNxy(2,:);
        B(3,1:2:7) = dNxy(2,:);
        B(3,2:2:8) = dNxy(1,:); 

        Ke = Ke + B'*C*B*wpg(ipg)*det(jac);
    end
    
    GdLEle = reshape(GdLNodos(Elementos(iele,:),:)',[],1);
    K(GdLEle,GdLEle) = K(GdLEle,GdLEle) + Ke;  
end


%% Rigid Links

CM = zeros(7,GdLTot);

u30 = GdLNodos(30,1);
u31 = GdLNodos(31,1);
u32 = GdLNodos(32,1);
u38 = GdLNodos(38,1);
u44 = GdLNodos(44,1);

v30 = GdLNodos(30,2);
v31 = GdLNodos(31,2);
v32 = GdLNodos(32,2);
v38 = GdLNodos(38,2);
v44 = GdLNodos(44,2);

CM(1,[u30 u38 u44])=[1 -2 1];
CM(2,[v31 v30 u30 u38])=[1 -1 -1 1];
CM(3,[v38 v30])=[1 -1];
CM(4,[v44 v30])=[1 -1];
CM(5,[u31 u30])=[1 -1];
CM(6,[u32 u30])=[1 -1];
CM(7,[v32 v30 u30 u38])=[1 -1 -2 2];
% CM(8,[u45 u38 u30])=[1 -2 1];
% CM(9,[v45 v30 u38 u30])=[1 -1 2 -2];

Q = zeros(7,1);

constrained = false(GdLTot,1);
constrained([v31 v32 v38 v44 u31 u32 u44 ]) = true;
released = ~constrained;

nconst = length(find(constrained));
Kc = [ K CM' ; CM zeros(nconst) ];
R = [ R ; zeros(nconst,1) + Q];


%% Resolucion

% Reduccion
Fijo = reshape(bc',[],1);
Libre = ~Fijo;
Libre = [ Libre ; true(nconst,1) ];
cond( Kc(Libre,Libre) );        %% Numero de condicionamiento: Infinito = No Inversible (o sea MAL)

% Solver
Dr = Kc(Libre,Libre)\R(Libre);

% Reconstruyo
D = zeros(GdLTot,1);
Libre = ~Fijo;
D(Libre) = D(Libre) + Dr(1 : end - nconst);

% Desplazamientos
D = (reshape(D,GdL,[]))';
CD = Nodos + D(:,1:2);


%% Graficación

figure('Name','Estructura Inicial y Desplazada','NumberTitle','off')
hold on
meshplot(Elementos,Nodos,'b')
meshplot(Elementos,CD,'r')
hold off
title('Estructura Inicial y Desplazada')


%% Salida de Datos

display('Las 9 fuerzas internas son:')
disp(' ')
disp(Dr(end-nconst+1 : end))

