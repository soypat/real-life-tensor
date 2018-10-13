tta=tita(2,:);
close all
%Potencia entregada a la masa
Potencias=figure;
    plot(tta,potenciaEntregadaALaMasa);
    hold on
    plot(tta,Pot)
    title('Potencia entregada al mecanismo y efectiva')
    ylabel('Potencia [W]')
    xlabel('\theta_2 [rad]')
    legend('Potencia Masa','Potencia Ent.')
    xlim([0 2*pi])
    
%Carga aplicada P
CargaAp=figure;
    plot(tta,Pv)
    title('Fuerzas aplicada sobre \alpha')
    ylabel('Fuerza [N]')
    xlabel('\theta_2 [rad]')
    xlim([0 2*pi])
    ylim([0 Pmax])
    text(3,25,'$$ F_{masa}=F_{max}(0,1+0,9 e^{-0,25|\theta-3|^{4} }) $$','Interpreter','latex')
%TORQUE
torque=figure;
    plot(tta,T)
    hold on
    plot(tta,Pot)
    title('Torque sobre manivela y potencia')
    ylabel('Torque/Potencia [Nm/W]')
    xlabel('\theta_2 [rad]')
    legend('Torque','Potencia')
    xlim([0 2*pi])
%FUERZAS
fuerzas=figure;
    plot(tta,A2)
    hold on
    plot(tta,A3)
    plot(tta,A6)
    plot(tta,B)
    plot(tta,C)
    plot(tta,D)
    plot(tta,O2)
    plot(tta,O4)
    plot(tta,Pv)
    legend('A2','A3','A6','B','C','D','O2','O4','P')
    title('Fuerzas absolutas sobre juntas')
    ylabel('Fuerza [N]')
    xlabel('\theta_2 [rad]')
    xlim([0 2*pi])
%ENERGIAS
potencias=figure;
    plot(tta,Ecinet(2,:));
    hold on
    plot(tta,Ecinet(3,:));
    plot(tta,Ecinet(4,:));
    plot(tta,Ecinet(5,:));
    plot(tta,Ecinet(6,:));
    legend('2','3','4','5','6')
    title('Energia cinetica de eslabones')
    xlabel('\theta_2 [rad]')
    ylabel('Energía neta [J]')
    