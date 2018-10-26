function Z=solucion_analitica2(X,Y)
%Devuelve la solución analítica del punto 2 evaluada en los puntos de las
%matrices X e Y. 

disp('Ni se te ocurra evaluarme afuera del dominio!')

L=0.5;%m
Tm=20;%C
b=1;%m  
    
Z=Tm*sinh(pi*Y/L)/sinh(pi*b/L).*sin(pi*X/L);