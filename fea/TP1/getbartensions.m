function [sig,N] = getbartensions(b,phi,Fv)
%GETVIGATENSIONS devuelve tensiones maximas para unaviga
s=sind(-phi);
c=cosd(-phi);
N=-(Fv(1)*c-Fv(2)*s);
A=b^2*pi/4;
sig=N/A;
end

