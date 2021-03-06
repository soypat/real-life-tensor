% clear all;
clc;
graph=false;
vueltas=3;

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
L5=L3;

BA=L3;
BD=161;
AX=186.5*escala;
N=120;
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
CA=norm(nodC(1,:)-nodA(1,:));
% hold on
refresh=true;

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
% h=gcf;

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

tita1=zeros(1,length(tita2));
angv1=zeros(length(angv2),1);
anga1=zeros(length(anga2),1);
anga=[anga1'; anga1'; anga3'; anga4']; %accel ang 2 es cero! Velocidad constante!
vang=[angv1';angv2';angv3';angv4'];
tita=[tita1;tita2;tita3;tita4];
%CAMBIO DE nombre!! tita2 = tita(2,:), accelAngular3= anga(3,:) etc. blabla
%% Calculo de Fuerzas
kari=150/sqrt((72+40)^2+(107+61)^2); %Escala Ari

titacg2=0; %medido desde O2
rcg2=sqrt(40^2 + 61^2)*kari;

titacg3=0; %Medido desde A
rcg3B=L3/2;

titacg4=0; %Medido desde O4
rcg4=(BD+L4)/2;

titacg5=titacg3;
rcg5= L3/2;

titacg6=27/180*pi;
rcg6C=norm([311 34]);


    N=length(tita2);
    Rcg2=zeros(N,2);
    Rcg3=zeros(N,2);
    Rcg4=zeros(N,2);
    Rcg5=zeros(N,2);
    Rcg6=zeros(N,2);
%Necesito acceleraciones de los centros de masa
for i=1:length(tita2)

    [RO2, RA2]=getR(tita2(i),titacg2,rcg2,L2);
    [RA3, RB3]=getR(tita3(i),titacg3,rcg3B,L3);
    [RO4, RD4]=getR(tita4(i),titacg4,rcg4,L4+BD);
    [RC5, RD5]=getR(tita3(i),titacg5,rcg5,L3);
    [RC6, RA6]=getR(tita4(i)+pi,titacg6,rcg6C,CA);
    
    Rcg2(i,:)=-RO2;
    Rcg3(i,:)=-RA3+nodA(i,:);
    Rcg4(i,:)=-RO4+nodO4(i,:);
    Rcg5(i,:)=-RC5+nodC(i,:);
    Rcg6(i,:)=-RC6+nodC(i,:);
end
G=zeros(6,N);
G(2,:)=Rcg2(:,1)+1i*Rcg2(:,2);
G(3,:)=Rcg3(:,1)+1i*Rcg3(:,2);
G(4,:)=Rcg4(:,1)+1i*Rcg4(:,2);
G(5,:)=Rcg5(:,1)+1i*Rcg5(:,2);
G(6,:)=Rcg6(:,1)+1i*Rcg6(:,2);

vG=loopdiff(G,tau);
aG=loopdiff(vG,tau);


%% Graph that stuff
%comment continue to graph
for j=1:vueltas
    y=tic;
    t=tic;
for i=1:N
    if ~graph
        break
    end
%     continue
    
    while refresh
        pause(.00000001)
    if toc(t)>(.7*tau)
        refresh=false;
    end
    end
    
        [RO2, RA2]=getR(tita2(i),titacg2,rcg2,L2);
    [RA3, RB3]=getR(tita3(i),titacg3,rcg3B,L3);
    [RO4, RD4]=getR(tita4(i),titacg4,rcg4,L4+BD);
    [RC5, RD5]=getR(tita3(i),titacg5,rcg5,L3);
    [RC6, RA6]=getR(tita4(i)+pi,titacg6,rcg6C,CA);
    RB4=RD4/2.28872168; %Alargo el vector para que llegue a B
%         aria=[0 19738 ]
    RP6=RA6+AX*(RA6-RC6)/norm(RA6-RC6);
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
    hold on
    scatter(Rcg2(i,1),Rcg2(i,2),'+')
    scatter(Rcg3(i,1),Rcg3(i,2),'+')
    scatter(Rcg4(i,1),Rcg4(i,2),'+')
    scatter(Rcg6(i,1),Rcg6(i,2),'+')
    scatter(Rcg5(i,1),Rcg5(i,2),'+')
    %Verificacion posiciones de Rcg
%     scatter(real(G(2,i))+RO2(1),imag(G(2,i))+RO2(2),'x')
%     scatter(real(G(2,i))+RA2(1),imag(G(2,i))+RA2(2),'x')
    
%     scatter(real(G(3,i))+RA3(1),imag(G(3,i))+RA3(2),'x')
%     scatter(real(G(3,i))+RB3(1),imag(G(3,i))+RB3(2),'x')
%     [RO4, RB4]
%     [RC5, RD5]
%     [RC6, RA6]
%     scatter(real(G(4,i))+RO4(1),imag(G(4,i))+RO4(2),'x')
%     scatter(real(G(4,i))+RB4(1),imag(G(4,i))+RB4(2),'x')
%     scatter(real(G(4,i))+RD4(1),imag(G(4,i))+RD4(2),'x')

%     scatter(real(G(5,i))+RC5(1),imag(G(5,i))+RC5(2),'x')
%     scatter(real(G(5,i))+RD5(1),imag(G(5,i))+RD5(2),'x')
    
%     scatter(real(G(6,i))+RC6(1),imag(G(6,i))+RC6(2),'x')
%     scatter(real(G(6,i))+RA6(1),imag(G(6,i))+RA6(2),'x')
    scatter(real(G(6,i))+RP6(1),imag(G(6,i))+RP6(2),'x')
    

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
%     pause(1)
    refresh=true;
end
%     toc(t)
end


function [Rcgin, Rcgout] = getR(titaglob,titacg,rin,L)
    titatot=titaglob+titacg;
    Rcgin=-rin*[cos(titatot) sin(titatot)];
    Rcgout=Rcgin+L*[cos(titaglob) sin(titaglob)];
end


function [newvec] = loopdiff(vec,tau)
    [~, N] = size(vec);
    newvec=zeros(size(vec));
    for i=1:N
    if i==1
        dx=vec(:,1)-vec(:,N-1);

    else
        dx=vec(:,i)-vec(:,i-1);
    end
    newvec(:,i)=dx/tau;
    end
end