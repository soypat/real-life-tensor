function [sig,tau,N,M] = getvigatensions(b,h,Fv)
%GETVIGATENSIONS devuelve tensiones maximas para unaviga
% s=sind(-phi);
% c=cosd(-phi);
N=Fv(4);
V=Fv(2);
M1=abs(Fv(3));
M2=abs(Fv(6));
M=max(M1,M2);
A=b*h;
I=A*h^2/12;
sig=N/A+M*h/(2*I);
tau=3/2*V/A;
end

