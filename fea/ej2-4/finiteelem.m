function desp = finiteelem(longitud,A_i,A_s,mod_elas,N,aprox,carga)
%Funcion que devuelve los desplazamientos en funcion de la...
%longitud,area,nodos y carga
%longitud en mm.
%presion en MPa
%Fuerza en N.
if N==1
    error('Pone mas nodos');
else
    nod=0:(longitud/(N-1)):longitud;%coloco los nodos    
    elementos=zeros(N-1,2);
    area=(A_i*(longitud-nod)/longitud)+A_s*nod/longitud;%calculo las areas de los nodos
for i=1:N-1                         %ubico los nodos que ocupan mi barra
   elementos(i,1)=i;
   elementos(i,2)=i+1;
end
    mr_global=rig_global(N,mod_elas,aprox,area,elementos,longitud);
    bc=[0;ones(N-1,1)];
    R=[0;zeros(N-2,1); carga];
    mr_red=mr_global(find(bc),find(bc));
    mr_inv=mr_red^-1;
    desp=mr_inv*[zeros(N-2,1); carga];
end
end

function kglobal=rig_global(N,mod_elas,aprox,area,elementos,longitud)
k=[1 -1;-1 1];
kglobal=zeros(N,N);
    if strcmp(aprox,'Area maxima')
        for i=1:N-1
        klocal=((mod_elas*area(i))/(longitud/(N-1)))*k;
        kglobal(elementos(i,:),elementos(i,:))=kglobal(elementos(i,:),elementos(i,:))+klocal;
        end
    elseif strcmp(aprox,'Area minima')
        for i=1:N-1
        klocal=((mod_elas*area(i+1))/(longitud/(N-1)))*k;
        kglobal(elementos(i,:),elementos(i,:))=kglobal(elementos(i,:),elementos(i,:))+klocal;
        end
    else
        for i=1:N-1
        klocal=(mod_elas*(area(i+1)+area(i))*longitud/(2*(N-1)))*k;
        kglobal(elementos(i,:),elementos(i,:))=kglobal(elementos(i,:),elementos(i,:))+klocal;
        end
    end
end