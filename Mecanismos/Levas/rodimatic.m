function [paso_rodillo,circulo_de_base,circulo_primario] = rodimatic(s,radio_rodillo,radio_base_leva)
% Obtiene los puntos en espacio 2D donde pasa el centro del rodillo
N=length(s);
tita=0:2*pi/(N-1):2*pi;
Nt=length(tita);
paso_rodillo=zeros(Nt,2);
circulo_primario=paso_rodillo;
Rb=radio_base_leva+radio_rodillo;

for i=1:Nt
    ti=tita(i);
    circulo_primario(i,1)=cos(ti);circulo_primario(i,2)=sin(ti);
    rx=(Rb+s(i))*cos(ti);ry=(Rb+s(i))*sin(ti);
    paso_rodillo(i,1)=rx;paso_rodillo(i,2)=ry;
end

circulo_de_base=circulo_primario*radio_base_leva;
circulo_primario=circulo_primario*Rb;
end

