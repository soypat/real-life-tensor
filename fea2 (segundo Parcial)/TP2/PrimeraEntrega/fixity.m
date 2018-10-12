for i=1:size(nodos,1)
    if nodos(i,1)==0
        bc(i,:)=true;
    elseif nodos(i,2)==0
        bc(i,:)=true;
    elseif nodos(i,1)==110
        bc(i,:)=true;
    elseif nodos(i,2)==30 && (nodos(i,1)<=10||nodos(i,1)>=10)
        bc(i,:)=true;
    end
    % Fijo nodos del radio
%     nodoscirc=[10 60;100 60];
    dx1=nodos(i,1)-10;
    dx2=nodos(i,1)-100;
    dy=nodos(i,2)-60;
    L1=sqrt(dx1^2+dy^2);
    L2=sqrt(dx2^2+dy^2);
    if abs(L1-30)<1e-5 || abs(L2-30)<1e-5
        bc(i,:)=true;
    end
    if nodos(i,1)<0 || nodos(i,2)<0
        bc(i,:)=true;
    end
end