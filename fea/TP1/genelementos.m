function [elementos]=genelementos(elen,eletype,nod)
% GENELEMENTOS genera la matriz elementos a partir de una matriz de
% conexiones entre nodos.
    ndof=3;
    [N,~]=size(nod);
    hinges=length(eletype(eletype>2));
    Ndof=N*ndof;
    Ndoftot=Ndof+hinges;
    [Ne,~]=size(elen);
    elementos=zeros(Ne,6);
    vecnodo=@(n) [n*3-2 n*3-1 n*3];
    vecnodoh=@(n) [n*3-2 n*3-1];
    hingecount=0;%elementos = [1 2 3 4 5 6]
    for i=1:Ne
        switch eletype(i)
            case {1,2}
                elementos(i,:)=[vecnodo(elen(i,1)) vecnodo(elen(i,2))];
            case 3
                hingecount=hingecount+1;
                elementos(i,:)=[vecnodoh(elen(i,1)) Ndof+hingecount vecnodo(elen(i,1))];
            case 4
                hingecount=hingecount+1;
                elementos(i,:)=[vecnodo(elen(i,1)) vecnodoh(elen(i,1)) Ndof+hingecount];
        end
    end
end