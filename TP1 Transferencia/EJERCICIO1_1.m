 function []=EJERCICIO1_1(cantidad_de_simulaciones,precision_inicial,potencia_incrementos)
    for i=cantidad_de_simulaciones:-1:1
        EJERCICIO1_jr(precision_inicial*potencia_incrementos^(i-1));
    end
end


function [  ] = EJERCICIO1_jr( cantidad_de_volumenes )
%% Condiciones iniciales y de contorno
%Barra unidimensional, seccion circular,de longitud L, area A
%Tamb fija, temperatura a la izquierda=T0 fija, derecho aislado,
%% Considereaciones, convenciones y supocisiones
%Sistema de unidades SI 
%Propiedades independientes de la temperatura
%Estado estacionario
%Se hace uso de matrices esparsas
%% Condiciones
%Datos del problema
longitud=1; area_transversal=.1; k=1; calor_generado=10;
temperatura_izquierda=0; temperatura_ambiente=30; h_ambiente=10;
%Calculo de propiedades
diametro=(area_transversal*4/pi)^.5;
longitud_volumen=longitud/cantidad_de_volumenes;
area_longitudinal=pi*diametro*longitud_volumen;
volumen_finito=pi*longitud_volumen*diametro^2/4;
%% Resolucion por volumenes finitos
%Todos los volumenes se transfieren calor con los dos de al lado por
%conduccion y con ambiente por convenccion
diagonal=(-2*k*area_transversal/longitud_volumen-h_ambiente*area_longitudinal)*speye(cantidad_de_volumenes);
%El primer volumen se transfiere con la pared con temperatura fija, con de la
%derecha y por conveccion con el ambiente
diagonal(1,1)=(-3*k*area_transversal/longitud_volumen-h_ambiente*area_longitudinal);
%El ultimo volumen solo se transfiere con el elemento de la derecha y convecta al
%ambiente
diagonal(end,end)=(-k*area_transversal/longitud_volumen-h_ambiente*area_longitudinal);
%Todos los volumenes se transfieren calor por conduccion con el de la derecha
%(menos el primero)
diagonal_superior=[sparse(cantidad_de_volumenes-1,1),k*area_transversal/longitud_volumen*speye(cantidad_de_volumenes-1)];
diagonal_superior(end+1,:)=0;
%Todos los volumenes se transfieren calor por conduccion con el de la izquierda
%(menos el ultimo)
diagonal_inferior=[sparse(1,cantidad_de_volumenes-1);k*area_transversal/longitud_volumen*speye(cantidad_de_volumenes-1)];
diagonal_inferior(:,end+1)=0;
sistema_lineal=diagonal+diagonal_inferior+diagonal_superior;
%Todos los volumenes generan calor e intercambian calor con el ambiente por
%conveccion
condiciones=(-calor_generado*volumen_finito-h_ambiente*temperatura_ambiente*area_longitudinal)*(sparse(cantidad_de_volumenes,1)+1);
%El primer volumen tambien se intercambia con la pared
condiciones(1)=condiciones(1)-k*area_transversal/longitud_volumen*2*temperatura_izquierda;
campo_de_temperaturas=sistema_lineal^-1*condiciones;
%% Ploteamos el campo de temperaturas
h=figure('Name',['Analisis para  ',num2str(cantidad_de_volumenes),' volumenes'],'rend','painters','pos',[10 10 500 900]);
subplot(3,1,1)
plot(0:longitud_volumen:1-longitud_volumen,campo_de_temperaturas);
xlabel('Distancia del borde izquierdo (m)');
ylabel('Temperatura (Cº)');
title(['Perfil de temperatura en la barra para ',num2str(cantidad_de_volumenes),' volumenes']);
%% Calculamos el flujo de calor
flujo_de_calor_por_la_barra=diff(campo_de_temperaturas)*k;
flujo_de_calor_al_ambiente=h_ambiente*(campo_de_temperaturas-temperatura_ambiente);
%% Ploteamos el flujo de calor
subplot(3,1,2)
plot(0:longitud_volumen:1-2*longitud_volumen,flujo_de_calor_por_la_barra)
xlabel('Distancia del borde izquierdo (m)');
ylabel('Flujo de calor en la barra $\frac{W}{m^2}$','Interpreter','latex');
title(['Flujo de calor en la barra para ',num2str(cantidad_de_volumenes),' volumenes']);
subplot(3,1,3);
plot(0:longitud_volumen:1-longitud_volumen,flujo_de_calor_al_ambiente)
xlabel('Distancia del borde izquierdo (m)');
ylabel('Flujo de calor hacia el ambiente $\frac{W}{m^2}$','Interpreter','latex');
title(['Flujo de calor hacia el ambiente para ',num2str(cantidad_de_volumenes),' volumenes']);
end



