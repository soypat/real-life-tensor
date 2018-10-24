%% Voy a analizar la transferencia de calor sobre mi regla de acero inoxidable mota REG030
% 
clear
%% Dimensiones del problema
L=1; % metros
Area=.1; %m^2 es la sección
k0=1; %Watt per meter kelvin
k=k0;
D=sqrt(4*Area/pi); 
qvol=10; %W/m^3

%% Condiciones Borde
T0=0;
Tamb=30;
hc=10;

%% Begin problem
N=11;%mis nodos
Ne=N-1; %Elementos
Le=L/Ne; %long. elementos

Ti=zeros(N,1);
Ti(1)=T0;%Centigrado


nodos=0:L/Ne:L;
C=zeros(N);
Q=zeros(N,1);
k=k0;
qgen = qvol*Area;

a=get_a(i,N);

b=get_b_temp(1,1);
c=0;
d=get_d_temp1orN(T0);
C(1,1+1)=-b;

C(N,N-1)=-1;
C(N,N)=1;
    
for i=2:N-1
    % No consideramos calor almacenado (inercia termica) d/dt=0 R.E.
    b=get_b_conveccion(i,N,L,hc,k);
    c=get_c_conveccion(i,N,L,hc,k);
    d=get_d_conveccion(i,N,L,k,qgen);
    C(i,i-1)=-c;
    C(i,i+1)=-b;
    C(i,i)=a;
    Q(i)=d;
end
T=C\Q;


figure(1)
plot(nodos,T,nodos,0.0161*solucion_analitica1(nodos))
title('Perfil de Temperaturas')
ylabel('Temperatura')
xlabel('Posicion')


function [d] = get_d_conveccion(i,N,L,k,qvol)
    Le=L/(N-1);
    if i~=1 && i~=N
        d=Le^2/k*qvol;
    else
        d=NaN;
    end
end

function [c] = get_c_conveccion(i,N)
    if i~=1 && i~=N
        c=1;
    else
        c=NaN;
    end
end
function [d] = get_d_temp1orN(Ti)
    d=Ti;
end
function [b] = get_b_temp(i,N)
    if i~=1 && i~=N
        b=1;
        return
    else
        b=0;
    end
end

function [b] = get_b_conveccion(i,N,L,h,k)
    if i==1
        Le=L/(N-1);
        b=(1+h*Le/k)^-1;
    elseif i~=N
        b=1;
    else
        b=0;
        return
    end
end

function [a] = get_a(i,N)
    % GET_A devuelve A con (i,N)
    if i==1 || i==N
        a=1;
    else
        a=2;
        return
    end
end
