% function Du=linear_load(L, A, E, nele)
L=60; A=2; E=30*10^6; nele=3;
% Discretización
Tx=@(x) -10*x;          %Función de la carga
long=L/nele;            %Longitud de cada elemnto
nod=(0:long:L).';       %Posicion de cada nodo
nnod=length(nod);       %Número de nodos
% Matriz global
K=sparse(nnod,nnod);	%Matriz vacia
Ke=E*A/long*[1 -1; -1 1]; %Matriz local
for e=1:nele
    K([e,e+1], [e,e+1])=K([e,e+1], [e,e+1])+Ke;
end
% Criterio de carga
delt=Tx(long)*long/2;           %Triangulo de carga
inc=[1:2:2*nele-1];             %Vector incremental
area=delt*inc;                  %Vector de cargas en cada elemento
cargaizq=[area/2,0];            %Vector de carga en el nodo izq de cada elemento
cargader=[0,area/2];            %Vector de carga en el nodo der de cada elemento
sig=cargaizq+cargader;
sigr=(sig(1:end-1)).'; %Se elimina el ultimo nodo ya que es el que no tiene despazamiento por estar empotrado
Kr=K(1:end-1,1:end-1);
D=Kr\sigr;
Du=D(1)



% delt=Tx(long)*long/2;         %Triangulo de carga
% f=delt/2*ones(1,nnod);        %Triangulo por elemento a la mitad
% f2=f; f2(1)=0; f2(end)=0;
% g=delt*(0:1:nnod-1); g(end)=0;
% g2=delt*[0,0:1:nnod-2];        %Cuadrados por elemento
% sig=f+f2+g+g2;                %Carga total por elemento
% sig=sparse(1,nnod);
% for e=1:nele
%     p=(Tx(nod(e))+Tx(nod(e+1)))/2; %La carga media de cada elemento se divide en 2 y va a cada nodo del mismo
%     sig(e)=sig(e)+p/2;
%     sig(e+1)=sig(e+1)+p/2;
% end
% sig(end)=0;