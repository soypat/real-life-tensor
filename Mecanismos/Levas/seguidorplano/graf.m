% figure(2)
% plot(x,Rho_p);
% xlim([0 360])
% title('Radio de curvatura de paso \rho_p')
% ylabel('Radio de curvatura [mm]','fontsize',14)
% xlabel('Giro de leva[°]','fontsize',14)
% 
% figure(3)
% plot(x,T)
% xlim([0 360])
% title('Torque considerando roce')
% ylabel('Torque [Nmm]','fontsize',14)
% xlabel('Giro de leva[°]','fontsize',14)

% figure(4)
% subplot(2,1,1)
% plot(x,Fc)
% grid on
% xlim([0 360])
% title('Fuerza axial sobre seguidor')
% ylabel('Fc [N]','fontsize',14)
% xlabel('Giro de leva[°]','fontsize',14)
% subplot(2,1,2)
% plot(x,F)
% grid on
% xlim([0 360])
% title('Fuerza de contacto entre leva y rodillo')
% ylabel('F [N]','fontsize',14)
% xlabel('Giro de leva[°]','fontsize',14)

% figure(5)
% plot(x,Mvolt)
% xlim([0 360])
% grid on
% title('Momento de volteo')
% ylabel('Momento [Nm]','fontsize',14)
% xlabel('Giro de leva[°]','fontsize',14)

% figure(6)
% plot(x,phi)
% xlim([0 360])
% grid on
% title('\Phi para una vuelta completa')
% ylabel('\Phi [°]','fontsize',14)
% xlabel('Giro de leva[°]','fontsize',14)

%Gráfico desplazamiento
figure(1);
subplot(4,2,1:2)
plot(x,s);
title('Desplazamiento (s)','fontsize',14)
xlabel('Giro de leva[°]','fontsize',14)
ylabel('mm','fontsize',14)
grid on
grid minor
xlim([0 360])
ylim([-10 50])

%%Gráfico v

subplot(4,2,3)
plot(x,v)
title('Velocidad (v)','fontsize',14)
xlabel('Giro de leva[°]','fontsize',14)
ylabel('mm/rad','fontsize',14)
grid on
grid minor
xlim([0 360])

%%Gráfico a

subplot(4,2,5)
plot(x,a)
grid on
grid minor
title('Aceleracion (a)','fontsize',14)
xlabel('Giro de leva[°]','fontsize',14)
ylabel('mm/rad^2','fontsize',14)
xlim([0 360])

%%Gráfico j

subplot(4,2,7)
plot(x(vic1),j(vic1),'b-','LineWidth',1)
grid on
grid minor
hold on
plot(x(vic2),j(vic2),'b-','LineWidth',1)
% plot([b1 b1],[j1(b1+1) j2(1)],'r-')
% plot([b1+b2 b1+b2],[j2(b2+1) 0],'r-')
plot(x(vic3),j(vic3),'b-','LineWidth',1)
title('Jerk (j)','fontsize',14)
xlabel('Giro de leva[°]','fontsize',14)
ylabel('mm/rad^3','fontsize',14)
xlim([0 360])

%%Gráfico V
subplot(4,2,4)
plot(x,V);
grid on;
grid minor
title('Velocidad (V)','fontsize',14)
xlabel('Giro de leva[°]','fontsize',14)
ylabel('mm/s','fontsize',14)
xlim([0 360])

%%Gráfico A
subplot(4,2,6)
plot(x,A);
grid on
grid minor
title('Aceleracion (A)','fontsize',14)
xlabel('Giro de leva[°]','fontsize',14)
ylabel('mm/s^2','fontsize',14)
xlim([0 360])

%%Gráfico J
subplot(4,2,8)
plot(x(vic1),J(vic1),'b-','LineWidth',1)
hold on
plot(x(vic2),J(vic2),'b-','LineWidth',1)
grid on
grid minor
plot(x(vic3),J(vic3),'b-','LineWidth',1)
% plot([b1 b1],[J1(b1+1) J2(1)],'r-')
% plot([b1+b2 b1+b2],[J2(b2+1) 0],'r-')
title('Jerk (J)','fontsize',14)
xlabel('Giro de leva[°]','fontsize',14)
ylabel('mm/s^3','fontsize',14)
xlim([0 360])


vmax = max(v);
amax = max(a);
jmax = max(j);

vmaxpos = find(v==vmax,1,'first');
amaxpos = find(a==amax,1,'first');
jmaxpos = find(j==jmax,1,'first');

vmin = min(v);
amin = min(a);
jmin = min(j);

vminpos = find(v==vmin,1,'first');
aminpos = find(a==amin,1,'first');
jminpos = find(j==jmin,1,'first');

Vmax = vmax*w;
Amax = amax*(w^2);
Jmax = jmax*(w^3);

Vmin = vmin*w;
Amin = amin*(w^2);
Jmin = jmin*(w^3);


%% PHI













