function [] = verify(nod,elenod)
%   VERIFY   verifica uniones entre nodos.
[N,~]=size(nod);
[Ne,~]=size(elenod);
connectednodes=false(N,1);
for i = 1:Ne
    nod1=elenod(i,1);
    nod2=elenod(i,2);
    if nod1<1 || nod2<1 || nod1>N || nod2>N
        error('Error fatal verificando nodos.\n\n Verificar matriz elenod o matriz nod.')
    end
    connectednodes([nod1 nod2])=true;
end
if sum(connectednodes)~=N
    error('Nodo sin conectar. Verificar elenod.')
end


end
