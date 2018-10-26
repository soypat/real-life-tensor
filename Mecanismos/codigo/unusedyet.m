tita1=zeros(1,length(tita2));
angv1=zeros(length(angv2),1);
anga1=zeros(length(anga2),1);
anga=[anga1'; anga1'; anga3'; anga4']; %accel ang 2 es cero! Velocidad constante!
vang=[angv1';angv2';angv3';angv4'];
tita=[tita1;tita2;tita3;tita4];
%CAMBIO DE nombre!! tita2 = tita(2,:), accelAngular3= anga(3,:) etc. blabla
%% Calculo de Fuerzas
kari=150/sqrt((72+40)^2+(107+61)^2); %Escala Ari
ICG2 = 4580.141916; %mm^2
ICG3 = 66100.713840; %mm^2
ICG4 = 15890.378168;%mm^2
ICG5 = 33479.246276; %mm^2
ICG6 = 15890.378168; %mm^2

syms Ax Ay Bx By O2x O2y O4x O4y

titacg2=0; %medido desde O2
rcg2=sqrt(40^2 + 61^2)*kari;

titacg3=0; %Medido desde A
rcg3B=L3/2;

titacg4=0; %Medido desde B
rcg4=L4/2;

titacg6=27/180*pi;
rcg6C=sqrt(525^2+84^2)*kari;



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
    [RB4, RO4]=getR(tita4(i),titacg4,rcg4,L4);
    [RC6, RA6]=getR(tita4(i),titacg6,rcg6C,CA);
    Rcg2(i,:)=-RO2;
    
    Rcg3(i,:)=-RA3+nodA(i,:);
    Rcg4(i,:)=-RB4+nodB(i,:);
end
%La carga sobre la amasadora va ser P, con Px y Py
for i=1:length(tita2)
    %ECUACIONES PARA BARRAS
    %Barra2
    [RO2, RA2]=getR(tita2(i),titacg2,rcg2,L2);
    [RA3, RB3]=getR(tita3(i),titacg3,rcg3,L3);
    [RB4, RO4]=getR(tita4(i),titacg4,rcg4,L4);
        Px=2;
        m2=3;
        T=10;
        
%         Fx2 = Ax+O2x == m2*aCG2(1,i); %Tomo Fuerza en junta A como positiva sobre eslabon2
%         Fy2 = Ay+O2y == m2*aCG2(2,i);
%         M2  = ICG2*anga(2,i)==T+RO2(1)*O2y-RO2(2)*O2x+RA2(1)*Ay-RA2(2)*Ax;
%         %T= Torque sobre eslabon. P=omega*T
%         M3  = ICG3*anga(3,i)==RA3(1)*(-Ay)-RA3(2)*(-Ax)+RB3(1)*By-RB3(2)*Bx;
%         M4  = ICG4*anga(4,i)==RB4(1)*(-By)-RB4(2)*(-Bx)+ RO4(1)*O4y-RO4(2)*O4x
%         Fx3 = -Ax+Bx + Px == Px

    
end
% function [Rcgin, Rcgout] = getR(titaglob,titacg,rin,L)
%     titatot=titaglob+titacg;
%     Rcgin=-rin*[cos(titatot) sin(titatot)];
%     Rcgout=Rcgin+L*[cos(titaglob) sin(titaglob)];
%end


