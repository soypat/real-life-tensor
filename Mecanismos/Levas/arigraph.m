
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
