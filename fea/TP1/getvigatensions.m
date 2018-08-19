function [sig,tau] = getvigatensions(b,h,Fv)
%GETVIGATENSIONS devuelve tensiones maximas para unaviga
[N, V, M1, M2]=[Fv(1)*c-Fv(2)*s;c*Fv(2)+s*Fv(1);abs(Fv(3));abs(Fv(6))];
M=max(M1,M2);
A=b*h;
I=A*h^2/12;
sig=N/A+M*h/(2*I);
tau=3/2*V/A;
end

