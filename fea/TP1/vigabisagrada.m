function [hingota] = vigabisagrada(E,A,I,L,start)
% VIGABISAGRADA crea la matriz rigidez local para viga bisagrada
%   VIGABISAGRADA(E,A,I,L,1 o 0) 1=comienza con bisagra, 0 finaliza con
%   bisagra.
if start
    hingota=(3*E*I/L^3)*[A*L^2/(3*I) 0 0 -A*L^2/(3*I) 0 0;
    0 1 0 0 -1 L;
    0 0 0 0 0 0;
    -A*L^2/(3*I) 0 0 A*L^2/(3*I) 0 0;
    0 -1 0 0 1 -L;
    0 L 0 0 -L L^2];
    return

elseif start==0
    hingota=(3*E*I/L^3)*[A*L^2/(3*I) 0 0 -A*L^2/(3*I) 0 0;
    0 1 L 0 -1 0;
    0 L L^2 0 -L 0;
    -A*L^2/(3*I) 0 0 A*L^2/(3*I) 0 0;
    0 -1 -L 0 1 0;
    0 0 0 0 0 0];
    return
end
end

