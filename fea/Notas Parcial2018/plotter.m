clear
close all
load parcial1.mat
load parcial2.mat
N1=numel(parcial1);
N2=numel(parcial2);%Notas 2018 2Cuatri

Nt=max(N1,N2);

par1=[parcial1; nan(Nt-N1,1)];
par2=[parcial2; nan(Nt-N2,1)];%Otorgo ausencias como NaNs

par=[par1, par2];

boxplot(par)
title('Box-whisker plot de los parciales')