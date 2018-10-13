%% GRAFICOS
% figure
% subplot(2,1,1)
%% VELOCIRAPTOR
% plot(tita2,vabsA,'Color',colorA)
% hold on
% plot(tita2,vabsB,'Color',colorB)
% plot(tita2,vabsC,'Color',colorC)
% plot(tita2,vabsD,'Color',colorD)
% plot(tita2,vabsX,'Color',colorX)
% text(pi,100,'A','Color',colorA)
% text(pi,80,'B','Color',colorB)
% text(pi,60,'C','Color',colorC)
% text(pi,40,'D','Color',colorD)
% text(pi,20,'\alpha','Color',colorX)
% title('Velocidades Absolutas')
% ylabel('mm/s')
% xlabel('\theta_2')

%% AXELERATORS
% plot(tita2,aabsB,'--','Color',colorB)

%% CALCULO Velocidades angulares

%% COMIENZAN GRAFICOS Angulos
% plot(tita2,tita2,'Color',colorA)
% hold on
% plot(tita2,tita3,'Color',colorB)
% plot(tita2,tita4,'Color',colorC)
% text(.5,1,'\theta_2','Color',colorA)
% text(1,5.2,'\theta_3','Color',colorB)
% text(2,3.86,'\theta_4','Color',colorC)
% xlim([0 2*pi])
% title('Angulos vs. \theta_2')
% ylabel('rad')
% xlabel('\theta_2')

% figure
%% velocidades angulares GRAF
% plot(tita2,angv2,'Color',colorA)
% hold on
% plot(tita2,angv3,'Color',colorB)
% hold on
% plot(tita2,angv4,'Color',colorC)
% text(1,-1.75,'\omega_2','Color',colorA)
% text(2,-.9,'\omega_3','Color',colorB)
% text(4,-.8,'\omega_4','Color',colorC)
% xlim([0 2*pi])
% title('Velocidades Angulares')
% ylabel('rad/s')
% xlabel('$$\theta_2$$','Interpreter','latex')
%% Axelerators angulators
% figure
% plot(tita2,anga2,'Color',colorA)
% hold on
% plot(tita2,anga3,'Color',colorB)
% plot(tita2,anga4,'Color',colorC)
% text(1,-.2,'$$\dot{\omega}_2$$','Interpreter','latex')
% text(3.2,-.6,'$$\dot{\omega}_3$$','Interpreter','latex')
% text(pi,.9,'$$\dot{\omega}_4$$','Interpreter','latex')
% xlim([0 2*pi])
% title('Aceleraciones Angulares')
% ylabel('$$rad/s^2$$','Interpreter','latex')
% xlabel('$$\theta_2$$','Interpreter','latex')   

eta=vabsX./vabsA;
plot(tita2,eta);
hold on
plot([0 2*pi],[1 1])
xlim([0 2*pi]);

title('Relación de velocidades entre \alpha y A')
ylabel('$$\eta_{\alpha A}=\frac{|V_{\alpha}|}{|V_{A}|}$$','Interpreter','latex')
xlabel('$$\theta_2$$','Interpreter','latex')
