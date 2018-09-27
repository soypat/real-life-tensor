function meshplot(elementos,nodos,color)
% MESHPLOT  Graficador de mallas.
% 
% MESHPLOT(elementos,nodos,color)
% 
% Numeraci�n de los nodos de los elementos:
%  4---7---3
%  |       |
%  8   9   6
%  |       |
%  1---5---2
%  
% elementos: matriz de conectividades.
% nodos:     matriz de coordenadas nodales.
% color:     string que especifica el color para los bordes de los
%            elementos.

nNodos = size(elementos,2);
switch nNodos
    case {8,9}
        vNod = [1 5 2 6 3 7 4 8];
    otherwise
        vNod = 1:nNodos;
end

h1 = patch('Faces',elementos(:,vNod),'Vertices',nodos);
set(h1,'EdgeColor',color,'FaceColor','none');

set(gca,'XTick',[],'YTick',[],'XColor',[1 1 1],'YColor',[1 1 1])
daspect([1 1 1])

for n=1:size(nodos,1)
    text(nodos(n,1),nodos(n,2),num2str(n))
end

for e=1:size(elementos,1)
    text(sum(nodos(elementos(e,:),1))/4,sum(nodos(elementos(e,:),2))/4,num2str(e))
end