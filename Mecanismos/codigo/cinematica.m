% clear all;
clc;
colorO2=[100 70 0]/255;
colorO4=[100 80 40]/255;
colorX=[183 27 0]/255;
colorA=[0 110 255]/255;
colorB=[221 147 0]/255;
colorC=[0 160 150]/255;
colorD=[60 160 0]/255;
rpm=15;
vueltas=1;
hz=rpm/60;
radps2=2*pi*hz; %radianes por segundo
escala=3.5714;
L1 = 480.95; %n1 a n4
L2 = 150;    %n1 a n2
L3 = 326.64; %n2 a n3
L4 = 410.86; %n3 a n4
BD=161;
AX=186.5*escala;
N=50;
h=2*pi/N;
tita2=linspace(2*pi,0,N);
tau=h/radps2;  %me dice el tiempo caracteristico de la simulacion. Cuanto tiempo pasa entre cada iteracion


[tita3,~,tita4,~]=solve_four_bars(tita2,L2,L3,L4,L1);

n1=0;
n2=L2*exp(1i*tita2);
n3=n2+L3*exp(1i*tita3);
n4=L1;


nodO2=zeros(N,2);%Final
nodA=zeros(N,2);
nodB=zeros(N,2);
nodC=zeros(N,2);
nodD=zeros(N,2);
nodO4=zeros(N,2);
nodaux=zeros(N,2);
nodX=zeros(N,2);
nodO4(:,1)=L1;
%% Node position calculation

for i=1:N
    points=[n1,n2(i),n3(i),n4];
    %POSITIONS
    nodA(i,:)=[real(n2(i)) imag(n2(i))];
    nodB(i,:)=[real(n3(i)) imag(n3(i))];
    vec4=(nodB(i,:)-nodO4(i,:))/norm(nodB(i,:)-nodO4(i,:));
    nodD(i,:)=nodB(i,:)+BD*vec4;
    vec3=(nodA(i,:)-nodB(i,:))/norm(nodA(i,:)-nodB(i,:));
    nodC(i,:)=nodD(i,:)+L3*vec3;
    vec6=(nodA(i,:)-nodC(i,:))/norm(nodA(i,:)-nodC(i,:));
    nodX(i,:)=nodA(i,:)+AX*vec6;
    %GRAPH
%     continue %comentar si queres grafico
    
%     pbaspect(1.6*[4 1 1])%axis equal
end

hold on
refresh=true;
%% Graph that stuff
%comment continue to graph
for j=1:10
    y=tic;
    t=tic;
for i=1:N
%     continue
    
    while refresh
        pause(.00000001)
    if toc(t)>(.7*tau)
        refresh=false;
    end
    end
    cla
    plot(nodX(:,1),nodX(:,2));
%     text(600,0,sprintf('N=%0.0f',i))
    txt=sprintf(' \theta_2 =%0.2f',tita2(i));
    text(400,-20,txt)
    axis equal
    scatter([nodO2(i,1) nodO4(i,1)],[nodO2(i,2) nodO4(i,2)],'ko')
    hold on
    scatter(nodD(i,1),nodD(i,2));
    hold on
    scatter(nodX(i,1),nodX(i,2));
    hold on
    scatter(nodC(i,1),nodC(i,2))
    line([nodO2(i,1) nodA(i,1)],[nodO2(i,2) nodA(i,2)],'Color','k');%ESLABON 2
    line([nodA(i,1) nodB(i,1)],[nodA(i,2) nodB(i,2)],'Color','k'); %ESLABON 3
%     line([nodB(i,1) nodO4(i,1)],[nodB(i,2) nodO4(i,2)]); %ESLABON 4 (B-O4)
    line([nodD(i,1) nodO4(i,1)],[nodD(i,2) nodO4(i,2)],'Color','k'); %ESLABON 4 (O4D)
    line([nodC(i,1) nodD(i,1)],[nodC(i,2) nodD(i,2)],'Color','k'); %ESLABON 5
    line([nodC(i,1) nodA(i,1)],[nodC(i,2) nodA(i,2)],'Color','k');%ESLABON 6
    line([nodA(i,1) nodX(i,1)],[nodA(i,2) nodX(i,2)],'Color',[.5 .5 .5],'LineStyle','--');%ESLABON 6 hasta alpha
    xlim([-L2-250 L1+L2+150]);
    ylim([-L1-100 L2+450]);
    title('Representación de la amasadora animada')
    ylabel('Posición y [mm]');
    xlabel('Posición x [mm]')
    t=tic;
    refresh=true;
end
    toc(T)
end
% hold on
% comet(nodX(:,1),nodX(:,2))

%% find velocity
velocA=zeros(N,2);
velocB=zeros(N,2);
velocC=zeros(N,2);
velocD=zeros(N,2);
velocX=zeros(N,2);
for i=1:N
    if i==1
        dxA=nodA(i,1)-nodA(N-1,1);
        dyA=nodA(i,2)-nodA(N-1,2);
        dxB=nodB(i,1)-nodB(N-1,1);
        dyB=nodB(i,2)-nodB(N-1,2);
        dxC=nodC(i,1)-nodC(N-1,1);
        dyC=nodC(i,2)-nodC(N-1,2);
        dxD=nodD(i,1)-nodD(N-1,1);
        dyD=nodD(i,2)-nodD(N-1,2);
        dxX=nodX(i,1)-nodX(N-1,1);
        dyX=nodX(i,2)-nodX(N-1,2);
    else
        dxA=nodA(i,1)-nodA(i-1,1);
        dyA=nodA(i,2)-nodA(i-1,2);
        dxB=nodB(i,1)-nodB(i-1,1);
        dyB=nodB(i,2)-nodB(i-1,2);
        dxC=nodC(i,1)-nodC(i-1,1);
        dyC=nodC(i,2)-nodC(i-1,2);
        dxD=nodD(i,1)-nodD(i-1,1);
        dyD=nodD(i,2)-nodD(i-1,2);
        dxX=nodX(i,1)-nodX(i-1,1);
        dyX=nodX(i,2)-nodX(i-1,2);
    end
    velocA(i,:)=[dxA dyA]/tau;
    velocB(i,:)=[dxB dyB]/tau;
    velocC(i,:)=[dxC dyC]/tau;
    velocD(i,:)=[dxD dyD]/tau;
    velocX(i,:)=[dxX dyX]/tau;
end
%% acceleration
axelA=zeros(N,2);
axelB=zeros(N,2);
axelC=zeros(N,2);
axelD=zeros(N,2);
axelX=zeros(N,2);
for i=1:N
    if i==1
        dvxA=velocA(i,1)-velocA(N-1,1);
        dvyA=velocA(i,2)-velocA(N-1,2);
        dvxB=velocB(i,1)-velocB(N-1,1);
        dvyB=velocB(i,2)-velocB(N-1,2);
        dvxC=velocC(i,1)-velocC(N-1,1);
        dvyC=velocC(i,2)-velocC(N-1,2);
        dvxD=velocD(i,1)-velocD(N-1,1);
        dvyD=velocD(i,2)-velocD(N-1,2);
        dvxX=velocX(i,1)-velocX(N-1,1);
        dvyX=velocX(i,2)-velocX(N-1,2);
    else
        dvxA=velocA(i,1)-velocA(i-1,1);
        dvyA=velocA(i,2)-velocA(i-1,2);
        dvxB=velocB(i,1)-velocB(i-1,1);
        dvyB=velocB(i,2)-velocB(i-1,2);
        dvxC=velocC(i,1)-velocC(i-1,1);
        dvyC=velocC(i,2)-velocC(i-1,2);
        dvxD=velocD(i,1)-velocD(i-1,1);
        dvyD=velocD(i,2)-velocD(i-1,2);
        dvxX=velocX(i,1)-velocX(i-1,1);
        dvyX=velocX(i,2)-velocX(i-1,2);
    end
    axelA(i,:)=[dvxA dvyA]/tau;
    axelB(i,:)=[dvxB dvyB]/tau;
    axelC(i,:)=[dvxC dvyC]/tau;
    axelD(i,:)=[dvxD dvyD]/tau;
    axelX(i,:)=[dvxX dvyX]/tau;
end
h=gcf;

%% velocidades/accel absolutas
vabsA=(velocA(:,1).^2+velocA(:,2).^2).^.5;
aabsA=(axelA(:,1).^2+axelA(:,2).^2).^.5;
vabsB=(velocB(:,1).^2+velocB(:,2).^2).^.5;
aabsB=(axelB(:,1).^2+axelB(:,2).^2).^.5;
vabsC=(velocC(:,1).^2+velocC(:,2).^2).^.5;
aabsC=(axelC(:,1).^2+axelC(:,2).^2).^.5;
vabsD=(velocD(:,1).^2+velocD(:,2).^2).^.5;
aabsD=(axelD(:,1).^2+axelD(:,2).^2).^.5;
vabsX=(velocX(:,1).^2+velocX(:,2).^2).^.5;
aabsX=(axelX(:,1).^2+axelX(:,2).^2).^.5;

% figure
% subplot(2,1,1)
%% VELOCIRAPTOR
plot(tita2,vabsA,'Color',colorA)
hold on
plot(tita2,vabsB,'Color',colorB)
plot(tita2,vabsC,'Color',colorC)
plot(tita2,vabsD,'Color',colorD)
plot(tita2,vabsX,'Color',colorX)
text(pi,100,'A','Color',colorA)
text(pi,80,'B','Color',colorB)
text(pi,60,'C','Color',colorC)
text(pi,40,'D','Color',colorD)
text(pi,20,'\alpha','Color',colorX)
title('Velocidades Absolutas')
ylabel('mm/s')
xlabel('\theta_2')

%% AXELERATORS
% plot(tita2,aabsB,'--','Color',colorB)

%% CALCULO Velocidades angulares
angv2=zeros(N,1);
angv3=zeros(N,1);
angv4=zeros(N,1);
for i=1:N
    if i==1
        dtit2=tita2(i+1)-tita2(i);
        dtit3=tita3(i)-tita3(N-1);
        dtit4=tita4(i)-tita4(N-1);
    else
        dtit2=tita2(i)-tita2(i-1);
        dtit3=tita3(i)-tita3(i-1);
        dtit4=tita4(i)-tita4(i-1);
    end
    angv2(i)=dtit2/tau;
    angv3(i)=dtit3/tau;
    angv4(i)=dtit4/tau;
end
%% CALCULO Aceleraciones angulares

anga2=zeros(N,1);
anga3=zeros(N,1);
anga4=zeros(N,1);
for i=1:N
    if i==1
        dw2=angv2(i)-angv2(N-1);
        dw3=angv3(i)-angv3(N-1);
        dw4=angv4(i)-angv4(N-1);
    else
        dw2=angv2(i)-angv2(i-1);
        dw3=angv3(i)-angv3(i-1);
        dw4=angv4(i)-angv4(i-1);
    end
    anga2(i)=dw2/tau;
    anga3(i)=dw3/tau;
    anga4(i)=dw4/tau;
end

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

% eta=vabsA./vabsX;
% plot(tita2,eta);
% plot([0 2*pi],[1 1])
% xlim([0 2*pi]);
% 
% title('Relación de velocidades entre A y \alpha')
% ylabel('$$\eta_{A \alpha}=\frac{|V_{A}|}{|V_{\alpha}|}$$','Interpreter','latex')
% xlabel('$$\theta_2$$','Interpreter','latex')

