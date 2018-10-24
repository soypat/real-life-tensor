function MeshplotTrigMec(elementos,nodos,bc,color,verbatim)
% MESHPLOT  Graficador de mallas.
%
% MESHPLOT(elementos,nodos,color)
%
% Numeración de los nodos de los elementos:
%  3
%  | \    
%  6   5  
%  |     \
%  1---4---2
%
% elementos: matriz de conectividades.
% nodos:     matriz de coordenadas nodales.
% color:     string que especifica el color para los bordes de los
%            elementos.
vNod = [1:3];

h1 = patch('Faces',elementos(:,vNod),'Vertices',nodos);
set(h1,'EdgeColor',color,'FaceColor','none');

set(gca,'XTick',[],'YTick',[],'XColor',[1 1 1],'YColor',[1 1 1])
daspect([1 1 1])
hold on
if verbatim
    for n=1:size(nodos,1)
        if bc(n,1)==1
            fill([nodos(n,1);nodos(n,1)-0.06;nodos(n,1)-0.06],[nodos(n,2); nodos(n,2)-0.03; nodos(n,2)+0.03],[1;1;1])
        end
        if bc(n,2)==1
            if bc(n,2)==1
                fill([nodos(n,1);nodos(n,1)-0.03;nodos(n,1)+0.03],[nodos(n,2); nodos(n,2)-0.06; nodos(n,2)-0.06],[1;1;1])
            end
        end
        text(nodos(n,1),nodos(n,2),num2str(n),'FontSize',18,'Color','r')
    end
    
    for e=1:size(elementos,1)
        text(mean(nodos(elementos(e,:),1)),mean(nodos(elementos(e,:),2)),num2str(e),'FontSize',18,'Color','b')
    end
end