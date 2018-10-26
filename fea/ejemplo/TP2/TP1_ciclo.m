%% Programa por correr
%Para correr este programa, se necesita el solver C2_5.m
%En el grupo se hicieron 2 solvers y una verificaci?n con ADINA: El informe
%es conjunto pero se puede verificar un programa de cada uno, subido a la
%carpeta de campus correspondiente.



clear all
close all 
clc

nTest = 300;

areaVector = zeros(nTest,1);
relativeMaxSigma = zeros(nTest,1);
sigmaBeamMax = zeros(nTest,1);
sigmaBarMax = zeros(nTest,1);
relativeNodeDisplacement = zeros(nTest,1); % Siempre comparo como se mueve el nodo inferior


for(e = 1:nTest)
    
    data = C2_5(e);
    sigmaMaxBeam = data.sigma.beam;
    sigmaBar = data.sigma.bar;
    barD = data.displacements.bars;
    beamD = data.displacements.beams;
    h = data.geometry.height;
    A = h^2;
    areaVector(e) = A;
    sigmaBeamMax(e) = max(abs(sigmaMaxBeam));
    sigmaBarMax(e) = max(abs(sigmaBar));
    relativeMaxSigma(e) = sigmaBeamMax(e)/sigmaBarMax(e);
    relativeNodeDisplacement(e) = norm(beamD(6,1:2))/norm(barD(6,:));
    
end

figure
plot(areaVector,sigmaBarMax./sigmaBeamMax,'r-')
title('$\sigma_{max}$ relativa vs Area','Interpreter','latex','FontSize',16)
xlabel('Area de la seccion $(mm^{2})$','Interpreter','latex','FontSize',15)
ylabel('$\frac{\sigma_{barra}}{\sigma_{viga}}$','Interpreter','latex','FontSize',22)
grid on

figure
plot(areaVector,relativeNodeDisplacement)
title('Comparacion del desplazamiento del nodo inferior','Interpreter','latex','FontSize',16)
xlabel('Area de la seccion $(mm^{2})$','Interpreter','latex','FontSize',15)
ylabel('$\frac{Viga}{Barra}$','Interpreter','latex','FontSize',18)
grid on
