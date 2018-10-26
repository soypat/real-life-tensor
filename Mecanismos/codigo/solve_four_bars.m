function [tita31,tita32,tita41,tita42]=solve_four_bars(tita2,a,b,c,d)
%Resuelve las posiciones de un mecanismo de cuatro barras. Devuelve tita3 y
%tita4. 31+41 o 32+42 componen la configuración cruzada o abierta del
%mecanismo(pag. 130 Norton)


K1=d/a;
K2=d/c;
K3=(a.^2-b.^2+c.^2+d.^2)/(2*a*c);
K4=d/b;
K5=(c.^2-d.^2-a.^2-b.^2)/(2*a*b);
%%
A = cos(tita2) - K1 - K2*cos(tita2) + K3;
B = -2*sin(tita2);
C = K1 - (K2+1)*cos(tita2) + K3;

D = cos(tita2) - K1 + K4*cos(tita2) + K5;
E = -2*sin(tita2);
F = K1 + (K4-1)*cos(tita2) + K5;

G = sqrt(B.^2 - 4*A.*C);
H = sqrt(E.^2 - 4*D.*F);

%%
if all(imag(G)==0)
tita41 = 2*atan2((-B+G),(2*A));
tita42 = 2*atan2((-B-G),(2*A));%Una no sirve
else
    disp('No grashof')
end

if all(imag(G)==0)
tita31 = 2*atan2((-E+H),(2*D));
tita32 = 2*atan2((-E-H),(2*D));%Una no sirve
else
    disp('No grashof')
end
end