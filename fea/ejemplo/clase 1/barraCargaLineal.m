%function[u]=barraCargaLineal(nElementos,loadFunction)
clear all
close all
clc


%%Datos geom?tricos y del material
nElementos=1000; %borrar
area = 2; % en in2
totalLength = 60;  % en in
E = 3E6;  % psi

%% Coordenadas, elementos y funcion de carga

x1 = 1:1:nElementos;
x2 = 2:1:nElementos+1;

dofPerNode = 1;

elementos = [x1' x2'];

nNodos = nElementos+1;

nodos = linspace(0,totalLength,nNodos);

Px =@(x) 600-10*x;%loadFunction; % Carga variable dependiente de x


%% Construccion de la matriz global


kGenerica =[1 -1;-1 1];
kGlobal = sparse(dofPerNode*nNodos,dofPerNode*nNodos);
loadVector = sparse(nNodos,1);


for(i = 1:nElementos)
    
    longitudElemento=nodos(i+1)-nodos(i);
    k = ((area*E)/longitudElemento);
    kLocal =k*kGenerica;
    
    kGlobal([i,i+1],[i,i+1]) = kGlobal([i,i+1],[i,i+1]) + kLocal;
    
    q = (Px(nodos(i))+ Px(nodos(i+1)))/2;  % Carga promedio sobre el elemento
    loadVector(i) = loadVector(i) + q/2;    % Divido la carga sobre el elemento entre sus nodos
    loadVector(i+1) = loadVector(i+1) + q/2;  
end


%% Vector de cargas
    
kReducida = kGlobal(2:end,2:end);
loadVectorReducido = loadVector(2:end);
D=sparse(nNodos,1);
D=kReducida\loadVectorReducido;
u=D(end)
a=nElementos/u

plot(nodos,loadVector)
    
    
    
 





