clear all;
clc;
b1 = 75;
b2 = 105;

%%Constantes del 1
c01 = 0;
c11 = 0;
c21 = 0;
c31 = 140;
c41 = -105;
%%Constantes del 2
c02 = 35;
c12 = 0;
c22 = -210;
c32 = 280;
c42 = -105;
%%
t1 = 0:1:b1;
t2 = 1:1:b2;

x1 = 0:1:b1;
x2 = b1:1:(b1+b2);
x3 = (b1+b2):1:360;

x = 0:1:361;

%%Desplazamiento
s1 = c01 + c11*(t1/b1).^1 + c21*(t1/b1).^2 + c31*(t1/b1).^3 + c41*(t1/b1).^4;
s2 = c02 + c12*(t2/b2).^1 + c22*(t2/b2).^2 + c32*(t2/b2).^3 + c42*(t2/b2).^4;
s3 = x3.*0;
%%Velocidad
v1 = c11 + 2*c21*(t1/b1).^1 + 3*c31*(t1/b1).^2 + 4*c41*(t1/b1).^3;
v2 = c12 + 2*c22*(t2/b2).^1 + 3*c32*(t2/b2).^2 + 4*c42*(t2/b2).^3;
v3 = x3.*0;
%%Aceleracion
a1 = 2*c21 + 6*c31*(t1/b1).^1 + 12*c41*(t1/b1).^2;
a2 = 2*c22 + 6*c32*(t2/b2).^1 + 12*c42*(t2/b2).^2;
a3 = x3.*0;
%%Jerk
j1 = 6*c31 + 24*c41*(t1/b1).^1;
j2 = 6*c32 + 24*c42*(t2/b2).^1;
j3 = x3.*0;

s = [s1 s2 s3];
v = [v1 v2 v3];
a = [a1 a2 a3];
j = [j1 j2 j3];

%%Calculo de magnitudes reales
w = 9*pi;
%%VELOCIDAD
V = v.*w;
%%ACELERACION
A = a.*(w^2);
%%JERK
J = j.*(w^3);



%%Gráfico desplazamiento
figure(1);
plot(x,s);
title('Desplazamiento (s)','fontsize',14)
xlabel('Grados[°]','fontsize',14)
ylabel('mm','fontsize',14)
grid on
xlim([0 360])

%%Gráfico v
figure(2);
subplot(2,1,1)
plot(x,v)
title('Velocidad (v)','fontsize',14)
xlabel('Grados[°]','fontsize',14)
ylabel('mm/rad','fontsize',14)
grid on
xlim([0 360])

%%Gráfico V
subplot(2,1,2)
plot(x,V);
grid on;
title('Velocidad (V)','fontsize',14)
xlabel('Grados[°]','fontsize',14)
ylabel('mm/s','fontsize',14)
xlim([0 360])

%%Gráfico a
figure(3);
subplot(2,1,1)
plot(x,a)
grid on
title('Aceleracion (a)','fontsize',14)
xlabel('Grados[°]','fontsize',14)
ylabel('mm/rad^2','fontsize',14)
xlim([0 360])

%%Gráfico A
subplot(2,1,2)
plot(x,A);
grid on
title('Aceleracion (A)','fontsize',14)
xlabel('Grados[°]','fontsize',14)
ylabel('mm/s^2','fontsize',14)
xlim([0 360])

%%Gráfico j
figure(4);
subplot(2,1,1)
plot(x,j)
grid on
title('Jerk (j)','fontsize',14)
xlabel('Grados[°]','fontsize',14)
ylabel('mm/rad^3','fontsize',14)
xlim([0 360])

%%Gráfico J
subplot(2,1,2)
plot(x,J);
grid on
title('Jerk (J)','fontsize',14)
xlabel('Grados[°]','fontsize',14)
ylabel('mm/s^3','fontsize',14)
xlim([0 360])

vmax = max(v);
amax = max(a);
jmax = max(j);

vmin = min(v);
amin = min(a);
jmin = min(j);



Vmax = vmax*w;
Amax = amax*(w^2);
Jmax = jmax*(w^3);

Vmin = vmin*w;
Amin = amin*(w^2);
Jmin = jmin*(w^3);















