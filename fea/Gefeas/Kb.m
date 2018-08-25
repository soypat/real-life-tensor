function [klocal] = Kb(E,A,L)
% TB crea matriz transformación
%   TB(phi)
k=E*A/L;
klocal=[k -k;
    -k k];

% T=[cosd(phi) 0;
%     sind(phi) 0;
%     0 cosd(phi);
%     0 sind(phi)];
% K=T*klocal*T';

% Kv=zeros(6);
% Kv([1 2 4 5],[1 2 4 5])=K;
end

