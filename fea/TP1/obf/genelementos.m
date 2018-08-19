function [elementos]=genelementos(elen)
% GENELEMENTOS genera la matriz elementos a partir de una matriz de
% conexiones entre nodos.

    [Ne,~]=size(elen);
    elementos=zeros(Ne,6);
    vecnodo=@(n) [n*3-2 n*3-1 n*3];
    for i=1:Ne
        elementos(i,:)=[vecnodo(elen(i,1)) vecnodo(elen(i,2))];
    end
end