%TP2
k=1;  %W/mK
L=0.5;%m
b=1;  %m

Nx=10;
Ny=20;

midpoint=@(A) (A(1:end-1)+A(2:end))/2;
x=linspace(0,L,Nx);
y=linspace(0,b,Ny);

