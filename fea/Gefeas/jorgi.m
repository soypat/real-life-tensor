h=1.8;


L1=2.8;
L2=2.1;
A1=L1*h;
A2=L2*h;
ro=1000;
g=9.81;
F1=ro*g*h/2*A1/2;
F2=ro*g*h/2*A2/2;

M = F1*L1/4+F2*L2/4;
V = (F1+F2)/sqrt(2);

%AISI 1020
Sy = 330e6; % Pa fluencia
n=1.7; % factor seguridad

syms t real % espesor
I = h*t^3/12;
Q = h*t^2/8;

tV=solve(Sy/n == V/(t*h),t );
fprintf('Ante corte simple: t=%0.1fmm\n',double(tV*1000))
tV=solve(Sy/n == V*Q/I/h,t );
fprintf('Ante corte jouravsky: t=%0.1fmm\n',double(tV*1000))

tFlexion = solve(Sy/n == M*h/(I*2),t );
fprintf('Ante flexion: t=%0.1fmm\n',double(tFlexion*1000))

