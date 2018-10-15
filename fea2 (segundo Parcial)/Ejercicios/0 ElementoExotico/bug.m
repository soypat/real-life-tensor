syms x y real
X = [1 x y x^2 x*y];
A = [1 -1 -1  1  1
     1  1 -1  1 -1
     1  1  1  1  1
     1 -1  1  1 -1
     1  0  0  0 0];  
shapefuns = X/A;

N(1,1:2:2*length(shapefuns))=shapefuns;
N(2,2:2:2*length(shapefuns))=shapefuns; %Tiene la forma de las funciones de forma encontradas en el cook pg 206, ecuacion (6.2-2). Despues veo si me sirven
dN(1,1:2:2*length(shapefuns))=diff(shapefuns,x);
u=1;
