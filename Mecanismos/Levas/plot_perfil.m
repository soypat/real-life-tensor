
paso_leva=levard;
plot(paso_rodillo(:,1),paso_rodillo(:,2));
hold on; plot(circulo_de_base(:,1),circulo_de_base(:,2));
plot(circulo_primario(:,1),circulo_primario(:,2));

scatter(0,0,'k','x');
plot(paso_leva(:,1),paso_leva(:,2),'k','LineWidth',2);%Perfil de leva
title('Perfil de curvas características')
legend('Recorrido de Rodillo','Circulo Base','Circulo Primario','Eje de giro','Perfil Leva')

daspect([1 1 1]);

