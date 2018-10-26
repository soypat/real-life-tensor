function [D]=regendesplazamientos(U,CB)
Ndof=length(CB);
CBneg=~CB;
D=zeros(Ndof,1);
Ucount=0;
%REgenera el vector U con 0 en desplazamientos restringidos
for i=1:length(CB)
    if CBneg(i)
        Ucount=Ucount+1;
        D(i)=U(Ucount);
    else
        continue
    end
end
end