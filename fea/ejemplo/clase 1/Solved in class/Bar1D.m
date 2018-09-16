% Modelo de barras 1D axial. Area variable
clear all
clc

% Dimensiones
Ab = 100;          %Area base mm2
At = 25;           %Area tope mm2
L = 4000;          %Alto columna mm
P = -1000;         %Carga N
E = 210000;        %N/mm2

Ex_disp = L*P/(E*(Ab - At)) * log(Ab/At);      %Solución exacta para el desplazamiento del extremo

A_eq = (P/Ex_disp)*(L/E);                      %Àrea equivalente para la solución exacta con 1 elemento

Nele = [2,2,4,8,16,32,64,128];

Displacements = zeros(2,length(Nele));
for n=1:length(Nele)                      %Numero de Elementos
  Disp =  Solve_bar(Ab,At,L,P,Nele(n),'P');
  Displacements(1:2,n) = [Nele(n);Disp];
  Disp =  Solve_bar(Ab,At,L,P,Nele(n),'m');
  Displacements(3,n) = Disp;
  Disp =  Solve_bar(Ab,At,L,P,Nele(n),'M');
  Displacements(4,n) = Disp;
end

figure(1)
hold on
plot(Displacements(1,:),(Displacements(2,:) - ones(1,length(Nele))*Ex_disp)/abs(Ex_disp),'k',...
    Displacements(1,:),(Displacements(3,:) - ones(1,length(Nele))*Ex_disp)/abs(Ex_disp),'b',...
    Displacements(1,:),(Displacements(4,:) - ones(1,length(Nele))*Ex_disp)/abs(Ex_disp),'r')
title('Convergencia')
xlabel('Número de elementos')
ylabel('Error relativo')
grid