%% DINAMICA
graph=true;
%{
Pregunta para pablo
Está bien modelar la fuerza como viscosa?
Cuanto le doy para la fuerza maxima? Tengo que usar los 11 kW?
Que onda las ecuaciones para calcular fuerzas. Esta bien 
descomponer en tres fuerzas sobre el bulon?
%}
% for i=2:length(tita2)-1

%% Dimensiones
%espesor e
e=10;%mm
meterpermm=1e-3;
density=7.9e3; %kg/m^3 ACERO
densitymm=density*(meterpermm)^3; %Kg/mm^3
Areas = [0 19738 79625 17256 51610 17256];%mm^2
masa=e*Areas*densitymm; %MASA CALCULADA

ICG2 = 4580.141916; %mm^2    IG*alfa= Nm = R*ma= m^2*kg*s^-2
ICG3 = 15890.378168; %mm^2
ICG4 = 33479.246276; %mm^2
ICG5 = 15890.378168;%mm^2
ICG6 = 66100.713840; %mm^2   mm^2*s^-2

IG=[0 ICG2 ICG3 ICG4 ICG5 ICG6];
IG=IG.*masa*1e-6; % kg m^2


%% Solver
syms A2x A2y A3x A3y A6x A6y Bx By Cx Cy Dx Dy O2x O2y O4x O4y To %17 incognitas
%FUERZA VISCOSA MODELADA
fP=@(t).1+ 0.9*exp(-1/4*abs(t-3).^4);%Fuerza en función de tita2
Pmax=500;
%Declarando variables
forza=[];
N=size(G,2);
anga=[anga;anga(3,:);anga(4,:)];
Pv=zeros(length(tita(2,:)),1);
for i=1:N
    %ECUACIONES PARA BARRAS
    %Barra2
    [RO2, RA2]=getR(tita2(i),titacg2,rcg2,L2);
    [RA3, RB3]=getR(tita3(i),titacg3,rcg3B,L3);
    [RO4, RD4]=getR(tita4(i),titacg4,rcg4,L4+BD);
    [RC5, RD5]=getR(tita3(i),titacg5,rcg5,L3);
    [RC6, RA6]=getR(tita4(i)+pi,titacg6,rcg6C,CA);
    normAX=AX;
    dirmov=velocX(i,:)/norm(velocX(i,:)); %Dirección de movimiento del punto de interes X
    P=Pmax*fP(tita(2,i))*(-dirmov); %Modelo de fuerza viscosa
    Pv(i)=norm(P);
    Px=P(1);Py=P(2);
    RP6=RA6+AX*(RA6-RC6)/norm(RA6-RC6);
    RB4=RD4/2.28872168;
        
    
        
%         Py=0;

%         T=1e3*15*746/abs(vang(2,i)); %Nmm

        %% NEWTON
        %Quiero fuerzas en Newtons: HAgo conversiones
%         IG=IG*1e-6;
        aG=aG/1000; %Nm %mm->m
        RO2=RO2/1000; RA2=RA2/1000; RA3=RA3/1000; RB3=RB3/1000; RO4=RO4/1000;
        RD4=RD4/1000; RC5=RC5/1000; RD5=RD5/1000; RC6=RC6/1000; RA6=RA6/1000; RB4=RB4/1000; RP6=RP6/1000; 
        %Defino B sobre eslabon 3, C sobre eslabon 6, D sobre eslabon 5
        Fx2 = A2x + O2x == masa(2)*real(aG(2,i)); %Tomo Fuerza en junta A como positiva sobre eslabon2
        Fy2 = A2y + O2y == masa(2)*imag(aG(2,i));
        M2  = To + RO2(1)*O2y-RO2(2)*O2x+RA2(1)*(A2y)-RA2(2)*A2x == 0;
        % Eslabon 3
        Fx3 = A3x + Bx  == masa(3)*real(aG(3,i));
        Fy3 = A3y + By  == masa(4)*imag(aG(3,i));
        M3  = RA3(1)*(A3y)-RA3(2)*(A3x)+RB3(1)*By-RB3(2)*Bx == IG(3)*anga(3,i);
        %Fuerzas sobre bulon (no tiene masa el buloncito. es de airgel.)
        FAx = A3x+A2x+A6x==0;
        FAy = A3y+A2y+A6y==0;
        %Eslabon 4
        Fx4 = -Bx + O4x - Dx == masa(4)*real(aG(4,i));
        Fy4 = -By + O4y - Dy == masa(4)*imag(aG(4,i));
        M4  =  RB4(1)*(-By)-RB4(2)*(-Bx)+ RO4(1)*O4y-RO4(2)*O4x+RD4(1)*(-Dy)-RD4(2)*(-Dy)==IG(4)*anga(4,i);
        %Eslabon 5
        Fx5 = -Cx + Dx == masa(5)*real(aG(5,i));
        Fy5 = Dy - Cy == masa(5)*imag(aG(5,i));
        M5  = RC5(1)*(-Cy)-RC5(2)*(-Cx) + RD5(1)*Dy - RD5(2)*Dx == IG(5)*anga(5,i);
        %T= Torque sobre eslabon. P=omega*T
        %eslabon 6 
        Fx6 = Px + A6x + Cx == masa(6)*real(aG(6,i));
        Fy6 = A6y + Cy + Py == masa(6)*imag(aG(6,i));
        M6  = - RP6(2)*Px +RP6(1)*Py + RC6(1)*Cy-RC6(2)*Cx + RA6(1)*A6y-RA6(2)*A6x == IG(6)*anga(6,i);
    %% Resuelvo NEWTON
%         M6  = RP6(1)*Py - RP6(2)*Px + RC6(1)*Cy-RC6(2)*Cx + RA6(1)*A6y-RA6(2)*A6x == ICG6*anga(6,i);
    equations=[Fx2 Fy2 M2 Fx3 Fy3 M3...
             Fx4 Fy4 M4 Fx5 Fy5 M5 Fx6 Fy6 M6 FAx FAy];   
    variables=[A2x A2y A3x A3y A6x A6y Bx...
            By Cx Cy Dx Dy O2x O2y O4x O4y To];
        
%         normalVar = [A2x A2y A3x A3y A6x A6y Bx By Cx Cy Dx Dy O2x O2y O4x O4y];
        [A,B] = equationsToMatrix(equations,variables);
    [X,R]=linsolve(A,B);

    Xe=eval(X);
    forza=[forza;Xe'];
    
end
%% Tomo Modulos
%variables=[A2x A2y A3x A3y A6x A6y Bx By Cx Cy Dx Dy O2x O2y O4x O4y To];
    A2=(forza(:,1).^2+forza(:,2).^2).^.5;
    A3=(forza(:,3).^2+forza(:,4).^2).^.5;
    A6=(forza(:,5).^2+forza(:,6).^2).^.5;
    B=(forza(:,7).^2+forza(:,8).^2).^.5;
    C=(forza(:,9).^2+forza(:,10).^2).^.5;
    D=(forza(:,11).^2+forza(:,12).^2).^.5;
    O2=(forza(:,13).^2+forza(:,14).^2).^.5;
    O4=(forza(:,15).^2+forza(:,16).^2).^.5;
    T=forza(:,17);
    
    %% Calculos de Energía
    Pot=T*vang(2,2);
    potenciaEntregadaALaMasa=abs(Pv.*vabsX/1000);
    eficienciaPot=potenciaEntregadaALaMasa./Pot;
    
    etaPot=potenciaEntregadaALaMasa./Pot;
    
    Ecinet=(abs(vG/1000).^2).*masa';
   %% Grafico
if graph==true

    dinagraf
end
%% OTRO




function [Rcgin, Rcgout] = getR(titaglob,titacg,rin,L)
    titatot=titaglob+titacg;
    Rcgin=-rin*[cos(titatot) sin(titatot)];
    Rcgout=Rcgin+L*[cos(titaglob) sin(titaglob)];
end


%Util para debugear
%     fprintf('Soluciónes:\n-----------------------\n')
%     for i=1:length(variables)
%         varstr=sprintf('%s',variables(i));
%         valstr=sprintf('%0.1f',Xe(i));
%         fprintf('%s = %s\n\n',varstr,valstr);
%     end
%     sol = solve([Fx2 Fy2 M2 Fx3 Fy3 M3...
%              Fx4 Fy4 M4 Fx5 Fy5 M5 Fx6 Fy6 M6 FAx FAy],...
%              [A2x A2y A3x A3y A6x A6y Bx By Cx Cy Dx Dy O2x O2y O4x O4y])
%     sol = solve([Fx2 Fy2 M2 Fx3 Fy3 M3...
%              Fx4 Fy4 M4 Fx5 Fy5 M5 Fx6 Fy6 M6 FAx FAy],...
%              [A2x A2y A3x A3y A6x A6y Bx By Cx Cy Dx Dy O2x O2y O4x O4y Px Py])
         
                  % end