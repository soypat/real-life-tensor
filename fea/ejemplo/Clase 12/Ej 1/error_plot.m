%% Error plot
clear; close all; clc
%% Selección de problema
eleT = 'Q4'; %'Q4', 'Q8', 'Q9', 'Q16', 'CST, 'LST'
nMallas = 5; % Cantidad de mallas

%% Carga de mallas
nodes1 = load('nodQ4_1.txt');
nodes2 = load('nodQ4_2.txt');
nodes3 = load('nodQ4_3.txt');
nodes4 = load('nodQ4_4.txt');
nodes5 = load('nodQ4_5.txt');
elements1 = load('eleQ4_1.txt');
elements2 = load('eleQ4_2.txt');
elements3 = load('eleQ4_3.txt');
elements4 = load('eleQ4_4.txt');
elements5 = load('eleQ4_5.txt');

%% Solución FEA
error = zeros(nMallas,2);

error(1,:) = [size(elements1,1) feasolve(eleT,nodes1,elements1)];
error(2,:) = [size(elements2,1) feasolve(eleT,nodes2,elements2)];
error(3,:) = [size(elements3,1) feasolve(eleT,nodes3,elements3)];
error(4,:) = [size(elements4,1) feasolve(eleT,nodes4,elements4)];
error(5,:) = [size(elements5,1) feasolve(eleT,nodes5,elements5)];

plot(error(:,1),error(:,2))