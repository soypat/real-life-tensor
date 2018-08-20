function [sig,tau,N,M] = getvigatensions(b,h,phi,Fv)
%GETVIGATENSIONS devuelve tensiones maximas para unaviga
s=sind(-phi);
c=cosd(-phi);
N=Fv(1)*c-Fv(2)*s;
V=c*Fv(2)+s*Fv(1);
M1=abs(Fv(3));
M2=abs(Fv(6));
M=max(M1,M2);
A=b*h;
I=A*h^2/12;
sig=N/A+M*h/(2*I);
tau=3/2*V/A;
end

