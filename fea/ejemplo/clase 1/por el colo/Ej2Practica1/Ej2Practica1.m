clear all
close all
clc

nMaxElementos = 1000;
elementos = 1:1:nMaxElementos;
u = sparse(nMaxElementos);
sigma = sparse(nMaxElementos);

for(i = 1:nMaxElementos)
    state = barraCargaLineal(i,linearLoad(0,-10));
    u(i) = state(1);
    sigma(i) = state(2);
    
end

figure(2)



plot(elementos,u,'b-')
xlabel('Numero de elementos')
ylabel('Desplazamiento')
grid


figure(3)
plot(elementos,sigma,'r*')
xlabel('Numero de elementos')
ylabel('Tension')


grid



    