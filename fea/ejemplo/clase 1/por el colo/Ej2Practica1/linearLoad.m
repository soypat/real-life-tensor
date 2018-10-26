% Toma un t?rmino independiente y una pendiente y saca una funci?n lineal
% dependiente de x con esos par?metros

function[P] = linearLoad(yIntercept,gradient)

P = @(x) yIntercept+gradient*x;
end