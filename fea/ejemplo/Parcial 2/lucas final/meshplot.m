function meshplot(elementos,nodos,color)
% MESHPLOT  Graficador de mallas.
% 
% MESHPLOT(elementos,nodos,color)
% 
% Numeración de los nodos de los elementos:
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

for i=1:size(elementos,1)
    verts = nodos(elementos(i,:),:);
    faces = [1 9 2 10 3 11 4 12; 5 17 6 18 7 19 8 20; 1 12 4 16 8 20 5 13; 2 10 3 15 7 18 6 14; 1 9 2 14 6 17 5 13; 3 11 4 16 8 19 7 15];
    patch('Faces',faces,'Vertices',verts,'FaceColor',color,'EdgeColor','w')
end
%  set(gca,'XTick',[],'YTick',[],'ZTick',[],'XColor',[1 1 1],'YColor',[1 1 1],'ZColor',[1 1 1])
daspect([1 1 1])
view(-30, 15);
end

