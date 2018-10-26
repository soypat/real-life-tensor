%% Ejercicio 1

% Para poder construir este elemento de 13 nodos sería necesario utilizar
% funciones de forma que contienen polinomios de orden 4. Se utilizarían 
% X = [ 1 x y x^2 x*y y^2 x^3 x^2*y x*y^2 y^3 x^3*y x^2*y^2 x*y^3] 
% excluyendo [x^4 y^4]. Esto significa que este elemento podrá
% representar polinomios cúbicos como máximo. Las tensiones y
% deformaciones representadas, al ser la derivada de dichos polinomios,
% por polinomios cuadráticos. Se le pueden aplicar cargas de orden cúbico.

% Mientras el jacobiano sea constante, la integral para conformar la matriz
% de rigidez tendra polinomios de orden 4, por lo tanto se necesita por lo
% minimo 3x3 puntos (2*3-1=5) para integrar full por regla de gauss.