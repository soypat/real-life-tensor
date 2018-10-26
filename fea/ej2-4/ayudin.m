
% DATOS (MEDIDAS EN mm)

E=210e3; I=0.5e8; A=0.5e4; L=[10e3 4e3]; q=300;

% DISCRETIZACION
DOF=3;
ELEMENTOS=[1 2;2 3;2 4;4 5];
NODOS=[0 0;0 4e3;10e3 4e3;0 8e3;10e3 8e3];
SIZE=size(ELEMENTOS);
Ne=SIZE(1);
Nn=Ne+1;
INDEX=ELEINDEX(DOF,ELEMENTOS);

% MATRIZ DE RIGIDEZ

KG=zeros(DOF*(Nn));
Y11=A*E/L(2)^3;Y12=12*E*I/(L(2)^2);Y13=4*E*I/(L(2));Y14=2*E*I/L(2);K1=A*E/L(2);
Y21=A*E/L(1)^3;Y22=12*E*I/(L(1)^2);Y23=4*E*I/(L(1));Y24=2*E*I/L(1);K2=A*E/L(1);
K1local=[Y11 0 -Y12 -Y11 0 -Y12;0 K1 0 0 -K1 0;-Y12 0 Y13 Y12 0 Y14;-Y11 0 Y12 Y11 0 Y12;...
    0 -K1 0 0 K1 0;-Y12 0 Y14 Y12 0 Y13];
K2local=[K2 0 0 -K2 0 0;0 Y21 Y22 0 -Y21 Y22;0 Y22 Y23 0 -Y22 Y24;-K2 0 0 K2 0 0;...
    0 -Y21 -Y22 0 Y21 -Y22;0 Y22 Y24 0 -Y22 Y23];
for i=1:Ne
    if mod(i,2)==0
        K=K2local
    else
        K=K1local
    end
   KG(INDEX(i,:),INDEX(i,:))= K + KG(INDEX(i,:),INDEX(i,:)) 
end

%Condiciones de borde

cb=~[1;1;1;0;0;0;1;1;1;0;0;0;1;1;1];
K_reducida=KG(cb,cb)
carga=[2820e3 0 -184e7 3e6 0 152e7 0 0 0 18e4 0 160e6 0 0 0]';
carga_reducida=carga(cb)
desp=K_reducida\carga_reducida