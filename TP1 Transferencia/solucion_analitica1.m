function vector_de_temperatura=solucion_analitica1(X)

%Evalua la solución analítica en el vector X
k=1;%W/(m.K)
L=1;%m
A=0.1;%m2
D=sqrt(4*A/pi);
qvol=10;%W/m3
T0=0;%C
Tamb=30;%C
hamb=10;%W/(m2.K)

syms T(x)
dT = diff(T,x);
d2T= diff(T,x,2);

cond1= T(0)==T0;
cond2= dT(L)==0;

ode = 0 == k*A*d2T + qvol*A + pi*D*hamb*(Tamb-T);

funcion_simbolica=dsolve(ode,[cond1,cond2]);

vector_de_temperatura=double(subs(funcion_simbolica,X));