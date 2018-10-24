%% Voy a analizar la transferencia de calor sobre mi regla de acero inoxidable mota REG030
% 
clear
%% Dimensiones del problema
L=0.32; % no tomo la parte agujereada
A=.5e-3*24e-3; %1mm x 24mm es la secci�n
k0=14.4;
alfa=0.387e-5;
rho=7817;
cp=461; %J/kg
%% Condiciones Borde
T0=100;
Tamb=25;

hc=10;
%% Begin problem

N=4; %mis nodos
Ne=N-1; %Elementos
Le=L/Ne; %long. elementos

Ti=zeros(N,1);
Ti(1)=T0;%Centigrado

nodos=0:L/Ne:L;
A=sparse(N);
Q=zeros(N,1);
k=k0;
qvol=2000;
for i=1:N
    x=nodos(i); 
    % No consideramos calor almacenado (inercia termica) d/dt=0 R.E.
    a=get_a(i,N);
    if i==1 %Solo en i==1 tenemos seteado temperatura
        b=get_b_temp(i,1);
        c=0;
        d=get_d_temp1orN(T0);
        A(i,i+1)=-b;
    elseif i~=N
        b=get_b_conveccion(i,N,L,hc,k);
        c=get_c_conveccion(i,N,L,hc,k);
        d=get_d_conveccion(i,N,L,hc,k,Tamb,qvol);
        A(i,i-1)=-c;
        A(i,i+1)=-b;
    else
        b=get_b_conveccion(i,N,L,hc,k);
        c=get_c_conveccion(i,N,L,hc,k);
        d=get_d_conveccion(i,N,L,hc,k,Tamb,qvol);
        A(i,i-1)=-c;
    end
    A(i,i)=a;
    Q(i)=d;
    
end

function [d] = get_d_conveccion(i,N,L,h,k,Tinf,qvol)
    Le=L/(N-1);
    if i~=1 || i~=N
        d=Le^2/k*qvol;
    else
        d=Le/k*(h*Tinf+qvol*Le/2)/(1+h*Le/k);
        return
    end
end

function [c] = get_c_conveccion(i,N,L,h,k)
    if i~=1 || i~=N
        c=1;
    elseif i==N
        Le=L/(N-1);
        c=(1+h*Le/k)^-1;
        return
    else
        c=0;
    end
end
function [d] = get_d_temp1orN(Ti)
    d=Ti;
end
function [b] = get_b_temp(i,N)
    if i~=1 || i~=N
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


   