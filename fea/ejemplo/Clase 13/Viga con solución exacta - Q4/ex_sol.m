function [ f ] = ex_sol(x,y,type)
%% Exacta
w=1;
P=80;
c=10;
L=100;
E = 1000;
NU = 0.25;
I=w*(2*c)^3/12;
G=E/(2+2*NU);

C = E/(1 - NU^2)*[ 1.0     NU         0.0
                    NU    1.0         0.0
                   0.0    0.0     (1 - NU)/2 ];
C_1 = inv(C);
if strcmp(type,'fun')
    f(:,1) = -(P*(x.^2-L^2).*y)/(2*E*I)-(NU*P.*y.*(y.^2-c^2))/(6*E*I)+(P.*y.*(y.^2-c^2))/(6*G*I);
    f(:,2) = (NU*P.*x.*y.^2)/(2*E*I)+(P*(x.^3-L^3))/(6*E*I)-((P*L^2)/(2*E*I)+(NU*P*c^2)/(6*E*I)+(P*c^2)/(3*G*I)).*(x-L);
elseif strcmp(type,'E')
    sx=(-3/2)*(P*x.*y)/c^3.*(C_1(1,1)*(-3/2)*(P*x.*y)/c^3);
    sy=0;
    txy=-(3*P/(4*c)).*(1-(y/c).^2).*C_1(3,3).* -(3*P/(4*c)).*(1-(y/c).^2);
    f = sx+sy+txy;
elseif strcmp(type,'stress')
    f = [(-3/2)*(P*x*y)/c^3;0;-(3*P/(4*c))*(1-(y/c)^2)];
end

