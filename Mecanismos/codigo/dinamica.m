%% SOLVER
syms A2x A2y A3x A3y A6x A6y Bx By Cx Cy Dx Dy O2x O2y O4x O4y To %16 incognitas
% for i=2:length(tita2)-1
fP=@(t).1+ 0.9*exp(-1/4*abs(t-3).^4);%Fuerza en función de tita2
Pmax=-200;

forza=[];
N=size(G,2);
anga=[anga;anga(3,:);anga(4,:)];
for i=1:N
    %ECUACIONES PARA BARRAS
    %Barra2
    [RO2, RA2]=getR(tita2(i),titacg2,rcg2,L2);
    [RA3, RB3]=getR(tita3(i),titacg3,rcg3B,L3);
    [RO4, RD4]=getR(tita4(i),titacg4,rcg4,L4+BD);
    [RC5, RD5]=getR(tita3(i),titacg5,rcg5,L3);
    [RC6, RA6]=getR(tita4(i)+pi,titacg6,rcg6C,CA);
    normAX=AX;
    dirmov=velocX(i,:)/norm(velocX(i,:));
    P=Pmax*fP(tita(2,i))*(-dirmov); %Modelo de fuerza viscosa
    Px=P(1);Py=P(2);
    RP6=RA6+AX*(RA6-RC6)/norm(RA6-RC6);
    RB4=RD4/2.28872168;
        
    
        
%         Py=0;
        m2=.3;
        m3=.3;
        m4=.3;
        m5=.3;
        m6=.3;
%         T=1e3*15*746/abs(vang(2,i)); %Nmm

        %NEWTON
        %Quiero fuerzas en Newtons: HAgo conversiones
        IG=IG*1e-6;
        aG=aG/1000;T=T/1000; %Nm %mm->m
        RO2=RO2/1000; RA2=RA2/1000; RA3=RA3/1000; RB3=RB3/1000; RO4=RO4/1000; RD4=RD4/1000; RC5=RC5/1000; RD5=RD5/1000; RC6=RC6/1000; RA6=RA6/1000; RB4=RB4/1000; RP6=RP6/1000; 
        %Defino B sobre eslabon 3, C sobre eslabon 6, D sobre eslabon 5
        Fx2 = A2x + O2x == m2*real(aG(2,i)); %Tomo Fuerza en junta A como positiva sobre eslabon2
        Fy2 = A2y + O2y == m2*imag(aG(2,i));
        M2  = To + RO2(1)*O2y-RO2(2)*O2x+RA2(1)*(A2y)-RA2(2)*A2x == 0;
        % Eslabon 3
        Fx3 = A3x + Bx  == m3*real(aG(3,i));
        Fy3 = A3y + By  == m3*imag(aG(3,i));
        M3  = RA3(1)*(A3y)-RA3(2)*(A3x)+RB3(1)*By-RB3(2)*Bx == IG(3)*anga(3,i);
        %Fuerzas sobre bulon (no tiene masa el buloncito. es de airgel.)
        FAx = A3x+A2x+A6x==0;
        FAy = A3y+A2y+A6y==0;
        %Eslabon 4
        Fx4 = -Bx + O4x - Dx == m4*real(aG(4,i));
        Fy4 = -By + O4y - Dy == m4*imag(aG(4,i));
        M4  =  RB4(1)*(-By)-RB4(2)*(-Bx)+ RO4(1)*O4y-RO4(2)*O4x+RD4(1)*(-Dy)-RD4(2)*(-Dy)==IG(4)*anga(4,i);
        %Eslabon 5
        Fx5 = -Cx + Dx == m5*real(aG(5,i));
        Fy5 = Dy - Cy == m5*imag(aG(5,i));
        M5  = RC5(1)*(-Cy)-RC5(2)*(-Cx) + RD5(1)*Dy - RD5(2)*Dx == IG(5)*anga(5,i);
        %T= Torque sobre eslabon. P=omega*T
        %eslabon 6 
        Fx6 = Px + A6x + Cx == m6*real(aG(6,i));
        Fy6 = A6y + Cy + Py == m6*imag(aG(6,i));
        M6  = - RP6(2)*Px +RP6(1)*Py + RC6(1)*Cy-RC6(2)*Cx + RA6(1)*A6y-RA6(2)*A6x == IG(6)*anga(6,i);
        
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

function [Rcgin, Rcgout] = getR(titaglob,titacg,rin,L)
    titatot=titaglob+titacg;
    Rcgin=-rin*[cos(titatot) sin(titatot)];
    Rcgout=Rcgin+L*[cos(titaglob) sin(titaglob)];
end