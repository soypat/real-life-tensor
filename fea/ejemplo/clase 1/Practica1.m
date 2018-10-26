%% Problema 1
% Nodos: Filas: #Nodos Columnas: Dim problema. Ej: 5 nodos en el plano ==> 5x2
% Elementos: Filas #Elementos Columnas: Nodos por elemento.
clear; close all; clc
L=4000; Af=100; Ai=25; E=210000; P=-1000;
A=@(x) 100-75*x/4000;
%% Estructura de datos geométricos
D=sparse(50,3);
for criterio=1:3 % Elección de área
for Nelem=1:50
Nodos=(0:L/Nelem:L).';
Nnodos=length(Nodos);
Elem=[1:Nelem;2:Nelem+1].';
areas=zeros(Nelem,1);
long=zeros(Nelem,1);
for e=1:Nelem
    if criterio==1
        areas(e)=A(Nodos(Elem(e,1)));
    else if criterio==2
            areas(e)=A(Nodos(Elem(e,2)));
        else 
            areas(e)=(A(Nodos(Elem(e,1)))+A(Nodos(Elem(e,2))))/2;
        end
    end
end
%% Matrices Locales
K=sparse(Nnodos,Nnodos);
for e=1:Nelem
    Ke=areas(e)*E/(L/Nelem)*[1 -1;-1 1];
    K([e,e+1],[e,e+1])=K([e,e+1],[e,e+1])+Ke;
end
%% Cargas
F=zeros(Nnodos,1);
F(Nnodos,1)=P;
Fr=F(2:end);
%% Condiciones de Borde
Kr=K(2:end,2:end);
%% Resolver
Dn=Kr\Fr;
D(Nelem,criterio)=Dn(end);
end
plot(1:50,D(:,criterio));
hold on
end
grid on



