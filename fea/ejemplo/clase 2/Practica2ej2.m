%% Solver de reticulados en 3D
clear; close all; clc
%% Datos de geométricos del problema
area=10000; 
young=210000;
F=180000;
%% Nodos
dofpornodo=3;       %Grados de libertad por nodo
nnod=8;             %Ingresar número de nodos
nod=1:1:nnod;       %Vector nombres de nodos
cnod=[-4500 0 3000
    4500 0 3000
    4500 0 -3000
    -4500 0 -3000
    -2250 10000 1500
    2250 10000 1500
    2250 10000 -1500
    -2250 10000 -1500]; %Coordenadas de cada nodo en filas
elem=[1 2
    1 3
    1 4
    1 5
    1 8
    2 3
    2 5
    2 6
    3 4
    3 6
    3 7
    4 7
    4 8
    5 6
    5 8
    6 7
    6 8
    7 8]; %Nodos vinculados por elementos representados en filas
nelem=size(elem,1); %Cantidad de elementos
%% Matriz de grados de libertad
dof=zeros(nnod,dofpornodo);
for n=1:nnod
    vec=(n-1)*dofpornodo+1:1:(n-1)*dofpornodo+3;
    dof(n,:)=vec;
end
%% Mariz global
Kg=zeros(dofpornodo*nnod);
A=area*ones(1,nelem); %Vector de areas
E=young*ones(1,nelem); %Vector de E
for e=1:nelem
    v=cnod(elem(e,2),:)-cnod(elem(e,1),:); %Vector de cada elemento
    long=norm(v);                          %Norma del vector
    vd=v/long;                             %Vector director
    T=[vd 0 0 0; 0 0 0 vd];                %Matriz de transformación
    Kl=A(e)*E(e)/long*[1 -1; -1 1];              %Matriz de rigidez local
    Ke=T.'*Kl*T;                           %Matriz de rigidez rotada del e.
    dofr=[dof(elem(e,1),:),dof(elem(e,2),:)]; %Vector de dof reducido al e.
    Kg(dofr,dofr)=Kg(dofr,dofr)+Ke;
end
%% Cargas y BC
Pr=[0 0 0 0 0 0 0 -F 0 F 0 0 F -F 0 0 -F 0].'; %Vector de cargas reducido (eliminando donde hay reacciones)
free=[4 6 9:1:24]; %Reductor de matriz global
Kgr=Kg(free,free);
Dr=Kgr\Pr;
%% Reacciones
D=[0 0 0 Dr(1) 0 Dr(2) 0 0 Dr(3:end).' ].';
P=Kg*D;
%% Tensiones
sig=zeros(nelem,1);
for e=1:nelem
    v=cnod(elem(e,2),:)-cnod(elem(e,1),:); %Vector de cada elemento
    long=norm(v);                          %Norma del vector
    vd=v/long;                             %Vector director
    elong=[D(elem(e,2)*3-2)-D(elem(e,1)*3-2) D(elem(e,2)*3-1)-D(elem(e,1)*3-1) D(elem(e,2)*3)-D(elem(e,1)*3)];
    sig(e)=E(e)/long*elong*vd.';
end
%% Máximo desplazamiento
for n=1:nnod
    Dn=D(n*3-2:1:n*3);
    Da(n)=norm(Dn);
end
pos1=find(abs(Da)==max(abs(Da)));
maxdesp=Da(pos1);
disp(strcat(['Nodo mas desplazado es el número: ',num2str(pos1)]))
disp(strcat(['Con un desplazamiento de: ',num2str(maxdesp),' mm']))
disp(' ')
%% Elemento mas tensionado
pos2=find(abs(sig)==max(abs(sig)));
maxtens=sig(pos2);
nod2=[elem(pos2,1) elem(pos2,2)];
disp(strcat(['Elemento mas tensionado es el número: ',num2str(pos2)]))
disp(strcat(['Entre los nodos: ',num2str(nod2(1)),' y ',num2str(nod2(2))]))
disp(strcat(['Con una tensión de: ',num2str(maxtens),' MPa']))

%% NOTA
% MUCHO CUIDADO CON EL REDUCTOR DE MA MATRIZ GLOBAL
    