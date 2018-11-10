clear all;
clc;
%{
%%Polinomio 1
A1 = [1,1;
     3,4]; 
B1 = [35;0];
C1 = linsolve(A1,B1);
%%Polinomio 2
A2 = [1,1,1;
     2,3,4;
     2,6,12];
B2 = [-35;0;0];
C2 = linsolve(A2,B2);
%}

%%Constantes del 1
c01 = 0;
c11 = 0;
c21 = 0;
c31 = 140;
c41 = -105;
%%Constantes del 2
c02 = 35;
c12 = 0;
c22 = -210;
c32 = 280;
c42 = -105;

b1 = 75;
b2 = 105;

t1 = 0:1:b1;
t2 = 1:1:b2;

x3 = (b1+b2):1:360;

%%Desplazamiento
s1 = c01 + c11*(t1/b1).^1 + c21*(t1/b1).^2 + c31*(t1/b1).^3 + c41*(t1/b1).^4;
s2 = c02 + c12*(t2/b2).^1 + c22*(t2/b2).^2 + c32*(t2/b2).^3 + c42*(t2/b2).^4;
s3 = x3.*0;
%%Velocidad
v1 = c11 + 2*c21*(t1/b1).^1 + 3*c31*(t1/b1).^2 + 4*c41*(t1/b1).^3;
v2 = c12 + 2*c22*(t2/b2).^1 + 3*c32*(t2/b2).^2 + 4*c42*(t2/b2).^3;
v3 = x3.*0;
%%Aceleracion
a1 = 2*c21 + 6*c31*(t1/b1).^1 + 12*c41*(t1/b1).^2;
a2 = 2*c22 + 6*c32*(t2/b2).^1 + 12*c42*(t2/b2).^2;
a3 = x3.*0;
%%Jerk
j1 = 6*c31 + 24*c41*(t1/b1).^1;
j2 = 6*c32 + 24*c42*(t2/b2).^1;
j3 = x3.*0;

s = [s1 s2 s3];
v = [v1 v2 v3];
a = [a1 a2 a3];
j = [j1 j2 j3];
vic1=1:length(s1);
vic2=(length(s1)+1):(length(s1)+length(s2));
vic3=(length(s1)+length(s2)+1):length(s);
%%Calculo de magnitudes reales
w = 9*pi;
%%VELOCIDAD
V = v.*w;
%%ACELERACION
A = a.*(w^2);
%%JERK
J = j.*(w^3);
J1 = j1.*(w^3);
J2 = j2.*(w^3);
J3 = j3.*(w^3);
x=linspace(0,360,length(s));
%Determinación de Rb

Rb = 0;

Rho_p = Rb+s+a;


while ~all(abs(Rho_p)>0)
   Rb = Rb + 1; 
   Rho_p = Rb+s+a;
end

% Rb = Rb + 2; %Haciendo esto el phi < 35°
[paso_rodillo,circulo_de_base,circulo_primario] = rodimatic(s,0,Rb);
levard = levamatic(paso_rodillo,0);%Estos valores para graficar
%Determinación de epsilon

cota = 0.0001;
phi = zeros(1,length(s));


% for e=-3:0.0001:3
%     if max(abs(phi))<phim
%         phim=max(abs(phi));
%         phig=phi;%Perfil hallado
%         em=e;
%     end
% end


% phi=phig;
n=0;
iter=-10:0.01:10;
Niter=length(iter);
M=zeros(Niter,1);
for e=iter
    n=n+1;
% fprintf('Fue encontrado un |Phi|=%0.2f \nPara una eccentricidad de %f\n',phim,em)

%El phimax es casi igual al phimin en módulo. Por lo tanto el epsilon es 0

%Determinacion del radio base para seguidor plano

amin = min(a);
aminpos = find(a==amin,1,'first');
smin = s(aminpos);

Rbplano = abs(amin + smin) + 1;
%% Fuerzas

newton_per_pound=0.453592*9.8;
inch_per_meter=(25.4E-3)^-1;

k=25*(newton_per_pound)*(inch_per_meter);
sigmoid=0.06; %Factor damping para levas generico
m=.45; %.5 kilogramos

x0=2.058; %mm 3.19
c=2*sigmoid*sqrt(m*k);
F=A*m+c*V+k*(s+x0);
Fc=F.*cosd(phi); %phi en grados

%% Si queres optimizar k, aca hay un calculin:
Fcminpos = find(Fc==min(Fc),1,'first');
Fcmin = Fc(Fcminpos);
sFcmin = s(Fcminpos);
k_tentativo = k-Fcmin/(s(Fcminpos)+x0);
x0_tentativo = x0-Fcmin/k;
Fcmax=max(Fc);
Fcmin/Fcmax %0.1

Mvolt=Fc.*(v/1000-e);
M(n)=max(abs(Mvolt));
end
plot(iter,M);
xlabel('Excentricidad [mm]')
ylabel('Momento de volteo [Nm]')
grid on
title('Momento de volteo maximo vs. excentricidad')
%% 
w_resonancia=sqrt(k/m-(c/(2*m))^2);
rpm=270;
verify

eta=1-.03;% 1 - Factor rozamiento
T=Fc.*V/(w*eta);

% plot(x,Fc)
%% Contacto curvo UNIDADES ESTANDAR
curv_leva=Rho_p/1000;
curv_rodillo=Rf/1000;
L=.2; %Longitud de leva metros
b=(  (2*F/(pi*L)).*( ( 1-.3^2 )/200e9+(1-.3^2)/200e9   )/(curv_rodillo^-1+curv_leva.^-1)  )^.5;
pmax=(2*F)./(b*pi*L);
% sigx=
