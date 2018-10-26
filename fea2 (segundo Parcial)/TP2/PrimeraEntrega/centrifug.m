    t=1;
    Diametro=280;
    R1=Diametro/2;%mm
    density=7.9e-3;%kg/cm^3
    density=density*(.1^3);%kg/mm^3
for iele=1:size(elementos,1)
    auxnodes=elementos(iele,:);
    auxpos=nodos(auxnodes,:);
    avgpos=sum(auxpos)/8;
    r=avgpos(2)+R1;
    accel=omega^2*r;
    m=density*Areas(iele)*t;
    Ftotal=m*accel/1000;%conversion de mN a N
    R(auxnodes(1:4),2)=R(auxnodes(1:4),2)-Ftotal/12*ones(4,1);
    R(auxnodes(5:8),2)=R(auxnodes(5:8),2)+Ftotal/3*ones(4,1);
end
%Busco nodos de pared sometida a compresion (y flexion)
paredzinha=[];
for i=1:nNod
    x=nodos(i,1);
    y=nodos(i,2);
    if (x==40)&& (y>=60)
        paredzinha=[paredzinha;y i];
    end
    
end

%HAgo sort de paredzinha
for i=1:size(paredzinha,1)
    for j=1:size(paredzinha)-i
        aux=paredzinha(j,:);
        aux2=paredzinha(j+1,:);
        if aux(1)<aux2(1)
            paredzinha(j+1,:)=aux;
            paredzinha(j,:)=aux2;
            continue
        end
    end
end
wallnodes=paredzinha(:,2);
% if paredzinha(1,1)>paredzinha(2,1)
%     aux=paredzinha(1,:);
%     paredzinha(1,:)=paredzinha(2,:);
%     paredzinha(2,:)=aux;
% end