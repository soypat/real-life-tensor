
% FUNCION PARA CONVERTIR ELEMENTOS EN VECTOR DE INDEXADO

function ELEDOF=ELEINDEX(DOF,ELEMENTOS)

% DOF=3;                            ejemplo...
% ELEMENTOS=[1 2;2 3;2 4;4 5];
SIZE=size(ELEMENTOS);
Ne=SIZE(1);
ELEDOF=zeros(Ne,DOF*2);

for i=1:Ne
    PRIMERNODO=ELEMENTOS(i,1)*DOF-(DOF-1):ELEMENTOS(i,1)*DOF;
    SEGUNDONODO=ELEMENTOS(i,2)*DOF-(DOF-1):ELEMENTOS(i,2)*DOF;
    ELEDOF(i,:)=[PRIMERNODO, SEGUNDONODO];
end
end