function [CB] = fixnode(CB,ndof,nodenumber)
%Empotra un nodo.
%fixnode(CondicionesBorde, grados de libertad por nodo, numero de nodo)
[m, n]=size(CB);
nodenumber=int64(nodenumber);
if n>1 || nodenumber*ndof>m
    error('lol bad CB')
end
if nodenumber<1
    error('bad node number')
end
vec=nodenumber:nodenumber+ndof;
CB(vec)=logical(vec);
return
end

