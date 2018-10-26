%% Rotador
function [Kp]=rotador(v,A,E);
long=norm(v);
K=A*E/long*[1 -1; -1 1];
b=v/long;
T=[b 0 0 0; 0 0 0 b];
Kp=T.'*K*T;