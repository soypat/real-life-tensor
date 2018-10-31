function [perfil_leva] = levamatic2(rodillo_path,Rr)
% toma los puntos del camino del rodillo y el Radio del rodillo
N=size(rodillo_path,1);
R=rodillo_path;
mdpts=zeros(N,2);
pleva=zeros(N,2);
%% Busco puntos medios
for i=1:N-1
    midp=(R(i,:)+R(i+1,:))/2;
    slope=(R(i+1,2)-R(i,2))/(R(i+1,1)-R(i,1));
    Rx=R(i+1,1)-R(i,1);Ry=R(i+1,2)-R(i,2);
    phi=atan2d(Ry,Rx);
    if slope == 0
        slope=eps;
    end
    phii=phi+pi/2;
    sector=sqrt((Rr/2)^2-(norm(R(i+1,:)-R(i,:)))^2);%distancia entre puntos es lo que esta adentro del norm()
    %Me tengo que fijar si está 'abajo' para sumarle en vez de restarle

    if ((R(i+1,1)-R(i,1))>0 && slope<0) || ((R(i+1,1)-R(i,1))>0 && slope>0)
        mdpts(i,:)=midp+[sector*cos(phii) sector*sin(phii)];
        
    else
        mdpts(i,:)=midp+[sector*cos(phii) sector*sin(phii)];
    end
end
%% Caso i==N
mdpts(N,:)=mdpts(1,:);
%% Busco perfil de leva con puntos medios
for i=1:N-1
    mdpt1=mdpts(i,:);mdpt2=mdpts(i+1,:);
    point=R(i+1,:);
    absmd=(mdpt1+mdpt2)/2;
    vec=(absmd-point)/norm(absmd-point);
    pleva(i,:)=Rr*vec+point;
end
%% Caso i==N
pleva(N,:)=pleva(1,:);
perfil_leva=pleva;
end

