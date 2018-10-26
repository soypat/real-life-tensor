function [B] = doper(N)
% operador mecanico
syms x y B real
col=size(N,2);
B=sym(zeros(3,col));
for i=1:3
    for j=1:col
        switch i
            case 1
                B(i,j)=diff(N(1,j),x);
            case 2
                B(i,j)=diff(N(2,j),y);
            case 3
                B(i,j)=diff(N(1,j),y)+diff(N(2,j),x);
        end
    end
end

end

