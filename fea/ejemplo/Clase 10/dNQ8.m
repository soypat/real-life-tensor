% FUNCION dNQ8
%Evalúa las derivadas de las funciones de forma de un Q8, respecto a csi y
%eta, en las coordenadas especificados en el input
% INPUT:    - csi, eta: coordenadas del punto dentro del elemento isoparamétrico
% OUTPUT:   - dN: matriz 2x8 con las derivadas de las funciones de forma

function dN = dNQ8(csi,eta)
dNcsi = [-((2*csi + eta)*(eta - 1))/4
    -((eta - 1)*(2*csi - eta))/4
    ((2*csi + eta)*(eta + 1))/4
    ((eta + 1)*(2*csi - eta))/4
    csi*(eta - 1)
    1/2 - eta^2/2
    -csi*(eta + 1)
    eta^2/2 - 1/2]';
dNeta = [ -((csi + 2*eta)*(csi - 1))/4
    -((csi - 2*eta)*(csi + 1))/4
    ((csi + 2*eta)*(csi + 1))/4
    ((csi - 2*eta)*(csi - 1))/4
    csi^2/2 - 1/2
    -eta*(csi + 1)
    1/2 - csi^2/2
    eta*(csi - 1)]';

dN = [dNcsi; dNeta];

