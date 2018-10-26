%% C2.13
clear all ; close all ; clc
%% Datos
Lbd = 1000; %mm
Lbc = 3000; %mm
Lab = 2000; %mm
Abd = 100; %mm^2
di = 20; %mm
de = 30; %mm
Ixy = pi()*(de^4-di^4)/64; %Momento de inercia para x e y de barras BC y AB mm^4
Ip = Ixy*2; %Momento polar (en este caso rigidez) para la barra AB mm^4
E = 210000; %N/mm^2
G = E/(2*(1+0.3));
%% Nodos y elementos
barras = 1; %Cantidad de barras
vigas = 10; %Cantidad de vigas
nelem = vigas+barras; %Número de elementos
cnodvigas = [zeros(1,vigas+1)
            0:Lbc/vigas:Lbc
            zeros(1,vigas+1)]; % Coordenadas de los nodos de las vigas
cnod = [0 0 Lbd
       cnodvigas.']; %Coordenadas de los nodos, origen en el nodo 2 (B)
nnod = size(cnod,1); %Número de nodos
elem = [1:1:nelem
        2:1:nnod]'; %Nodos que determinan cada elemento
%% DOF
% DOF para la barra
dofpornodobarra = 3;
dofbarra = (reshape([1:1:dofpornodobarra]',dofpornodobarra,1))';
% DOF para las vigas
dofpornodoviga = 3;
dofviga = (reshape([1:1:(nnod-1)*dofpornodoviga]'+dofbarra(end),dofpornodoviga,nnod-1))';
%% Matriz global
Kg = zeros(dofpornodobarra+dofpornodoviga*(nnod-1),dofpornodobarra+dofpornodoviga*(nnod-1));
% Soporte
v = cnod(elem(1,2),:)-cnod(elem(1,1),:);
vd = v/Lbd;
X = Abd*E/Lbd;
Kl = X*[1 -1; -1 1];
T = [vd 0 0 0; 0 0 0 vd];
Ke = T.'*Kl*T;
dofr = [dofbarra dofviga(1,1:3)];
Kg(dofr,dofr) = Kg(dofr,dofr)+Ke;
% Caño
for e = 2:nelem
    v = cnod(elem(e,2),:)-cnod(elem(e,1),:);
    long = norm(v);
    vd = v/long;
    Y4 = 2*E*Ixy/long;
    Y3 = Y4*2;
    Y2 = 3*Y4/long;
    Y1 = 2*Y2/long;
    S=  G*Ip/long;
    Kl = [Y1 Y2 0 -Y1 Y2 0
          Y2 Y3 0 -Y2 Y4 0
         -Y1 -Y2 S Y1 -Y2 -S
          Y2 Y4 0 -Y2 Y3 0
          0 0 0 0 0 0
          0 0 -S 0 0 S];
    Ke = rotar2D(Kl,cnod(elem(e,1),:),cnod(elem(e,2),:),[0 0 1]);
    dofr = [dofviga(elem(e-1,1),:),dofviga(elem(e-1,2),:)];
    Kg(dofr,dofr) = Kg(dofr,dofr)+Ke;
end
%% Cargas
P = zeros(nnod,10;
P(7) = 10000; %N
%% Reducción
fijos = zeros(1,dofpornodobarra+dofpornodoviga*(nnod-1));
fijos(1:3)=1;
fijos(end-2:end)