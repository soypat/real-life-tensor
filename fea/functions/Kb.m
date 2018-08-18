function [Kv] = Kb(E,A,L,phi)
% TB crea matriz transformaci�n
%   TB(phi)
k=E*A/L;
K=[k -k;
    -k k];

T=[cosd(phi) 0;
    sind(phi) 0;
    0 cosd(phi);
    0 sind(phi)];
K=T*K*T';

Kv=zeros(6);
Kv([1 2 4 5],[1 2 4 5])=K;
end

