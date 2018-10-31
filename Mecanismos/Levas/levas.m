clear all;
clc;
graph=false;
a = 5/12;

A = [a^3,a^4,a^5,a^6,a^7;               %s a 75°
     3*a^2,4*a^3,5*a^4,6*a^5,7*a^6;     %v a 75°
     1,1,1,1,1;                         %s a 180°
     3,4,5,6,7;                         %v a 180°
     6,12,20,30,42];                    %a a 180°
 
B = [35;0;0;0;0];                       %valores de (s,v,s,v,a respectivamente)
C = linsolve(A,B);                      %Constantes(c3,c4,c5,c6,c7)


x = 0:0.005:1;                           %x = tita/B1+B2

B1 = 5/12 * pi;
B2 = pi - B1;

tita = x.*(B1+B2);

s = C(1)*x.^3 + C(2)*x.^4 + C(3)*x.^5 + C(4)*x.^6 + C(5)*x.^7;
v = 3*C(1)*x.^3 + 4*C(2)*x.^4 + 5*C(3)*x.^5 + 6*C(4)*x.^6 + 7*C(5)*x.^7;
a1 = 6*C(1)*x.^3 + 12*C(2)*x.^4 + 20*C(3)*x.^5 + 30*C(4)*x.^6 + 42*C(5)*x.^7; 

if graph
figure(1);
plot(x,s);
figure(2);
plot(x,v);
figure(3);
plot(x,a1);
end
tita=0:2*pi/(length(s)-1):2*pi;
% tita=linspace(0,2*pi,length(s));
Rb=100;% Ni idea, puse un radio base cualquiera
Rr=28;%Radio rodillo
R=zeros(length(tita),2);
Rbase=zeros(length(tita),2);
for i=1:length(tita)
    ti=tita(i);
    Rbase(i,1)=Rb*cos(ti);Rbase(i,2)=Rb*sin(ti);
    rx=(Rb+s(i))*cos(ti);ry=(Rb+s(i))*sin(ti);
    R(i,1)=rx;R(i,2)=ry;
end
levard=levamatic(R,Rr); %Posiciones de la surpiperficie de la leva
%% Plotting
plot(R(:,1),R(:,2),'r');
hold on; plot(Rbase(:,1),Rbase(:,2),'g');
daspect([1 1 1]);
scatter(0,0,'k','x');
plot(levard(:,1),levard(:,2),'b');%Perfil de leva
legend('Recorrido de Rodillo','RadioBase Rodillo','Centro de eje','Perfil Leva')
