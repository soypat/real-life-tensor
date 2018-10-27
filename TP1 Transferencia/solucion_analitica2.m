function Z=solucion_analitica2(x,y)
%Devuelve la solución analítica del punto 2 evaluada en los puntos de las
%matrices X e Y. 
[X, Y] = meshgrid(x,y);
disp('Ni se te ocurra evaluarme afuera del dominio!')

L=0.5;%m
Tm=20;%C
b=1;%m  
    
Z=Tm*sinh(pi*Y/L)/sinh(pi*b/L).*sin(pi*X/L);

% Z=Tm*sinh(pi*Y/L).*sin(pi*X/L);

Z=flipud(Z);