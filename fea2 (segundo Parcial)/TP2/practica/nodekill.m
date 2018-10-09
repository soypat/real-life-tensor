function [newnodes, elementos] = nodekill(nodos,nodenumbering,elementos)
%NODEKILL Regenerates nodes
Nnode=size(nodos,1);
ndof=size(nodos,2);
nele=size(elementos,1);
elenodcount=size(elementos,2);
dictionary=[];
newnodes=[];
k=0;
dictionary=NaN(Nnode,1);
for i=1:Nnode
    isLegit=~isempty(find(elementos==i,1));
    if isLegit
       k=k+1;
       dictionary(nodenumbering(i))=k;
       newnodes=[newnodes;nodos(i,:)];
    end
end
for iele=1:nele
    for i=1:elenodcount
        elementos(iele,i)=dictionary(elementos(iele,i));
    end
end
end

