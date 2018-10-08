    t=60;
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