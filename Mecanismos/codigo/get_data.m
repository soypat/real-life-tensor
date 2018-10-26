function [out] = get_data(tita2,dato)
%Obtiene dato puntual para un tita2. 
%get_data(tita,vector_datos)
%Ejemplo get_data((180)*pi/180,velocD)
%get_data((180)*pi/180,vabsB)
    N=size(dato,1);
    interval=tita2/(2*pi);
    n=floor(N*interval);
    diff=N*interval-n;
    out=dato(N-n+1,:);
end

