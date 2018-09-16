function PlotStressCut(elements,Dnodes,avgStress8,Tension,Aspect,FilasEliminadas,ColumnasEliminadas,ElementosEspesor,ElementosH,ElementosLargo,Cond)

%% Elementos a no considerar
elem=[];
for i=1:ElementosH
    elem=[elem (1:FilasEliminadas)+ElementosLargo*(i-1)];
end
PrimerPiso=elem;
for j=2:ElementosEspesor
    elem=[elem PrimerPiso+ElementosH*ElementosLargo*(j-1)];
end

elem2=[];
for i=1:ColumnasEliminadas
    elem2=[elem2 (1:ElementosLargo)+ElementosLargo*(i-1)];
end
PrimerPiso=elem2;
for j=2:ElementosEspesor
    elem2=[elem2 PrimerPiso+ElementosH*ElementosLargo*(j-1)];
end

elementsOut=elements(elem,:);
elem=unique([elem elem2]);
elements(elem,:)=[];
avgStress8(elem,:,:)=[];
%% Plots




figure
title('Tensión en 8 Puntos de Gauss y promediada')
 bandplot(elements,Dnodes,avgStress8(:,:,Tension),[],[],20) 
 hold on
axis([0 32 -1 8 -1 3])
daspect([Aspect 1 1]);
view([0 0]);
xlabel('x')
ylabel('y')
zlabel('z')

% 
% figure
% title('Tensión en Extrapolada de 27 Puntos de Gauss y Promediada')
% bandplot(elements,Dnodes,AvgStress(:,:,Tension),[],[],20) % x
%  hold on
% if Cond
%  meshplot(elementsOut,Dnodes,[0.9 0.9 0.9]) 
% end
% daspect([Aspect 1 1]);
% view([0 0]);
% xlabel('x')
% ylabel('y')
% zlabel('z')
% 
% figure
% title('Tensión en Extrapolada de 27 Puntos de Gauss y Promediada con Tensión en Centros de Celda')
% bandplot(elements,Dnodes,AvgStressCentro(:,:,Tension),[],[],20) % x
%  hold on
% if Cond
%  meshplot(elementsOut,Dnodes,[0.9 0.9 0.9]) 
% end
% daspect([Aspect 1 1]);
% view([0 0]);
% xlabel('x')
% ylabel('y')
% zlabel('z')
end