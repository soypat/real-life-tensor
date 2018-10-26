function [K] = Kvuw(E,A,Iz,Iy,K,nu,L)
%KVUW Te hace la matriz gigante.
%   Kvuw(E,A,Iz,Iy,K,poisson,L)
    X=A*E/L;
    G=E/(2+2*nu);
    S=G*K/L;
    Y1=12*E*Iz/L^3;
    Y2=6*E*Iz/L^2;
    Y3=4*E*Iz/L;
    Y4=Y3/2;
    Z1=12*E*Iy/L^3;
    Z2=6*E*Iy/L^2;
    Z3=4*E*Iy/L;
    Z4=Z3/2;
    K=[.5*X 0 0 0 0 0,-X 0 0 0 0 0;
        0 .5*Y1 0 0 0 Y2,0 -Y1 0 0 0 Y2;
        0 0 .5*Z1 0 -Z2 0, 0 0 -Z1 0 -Z2 0;
        0 0 0 .5*S 0 0,0 0 0 -S 0 0;
        0 0 0 0 .5*Z3 0, 0 0 Z2 0 Z4 0;
        0 0 0 0 0 .5*Y3, 0 -Y2 0 0 0 Y4;
        ...
        0 0 0 0 0 0, .5*X 0 0 0 0 0;
        0 0 0 0 0 0, 0 .5*Y1 0 0 0 -Y2;
        0 0 0 0 0 0, 0 0 .5*Z1 0 Z2 0;
        0 0 0 0 0 0, 0 0 0 .5*S 0 0;
        0 0 0 0 0 0, 0 0 0 0 .5*Z3 0;
        0 0 0 0 0 0, 0 0 0 0 0 .5*Y3];
    K=K+K';%simetriza las cosas
    return
end