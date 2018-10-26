function [St] = Calcutodo(varargin)
%CALCUTODO claramente calcula todo
%ahora calculo las tensiones en un elemento viga
%tiene que ingresarle:
%1-todas las A.2-todas la I.3-todas las y.4-todas la b.5-todas las fl
%St=Calcutodo(A,I,y,b,E,fl1);

n=nargin;

A=[varargin{1:n/6}];
I=[varargin{(n/6+1:(2*n)/6)}];
y=[varargin{((2*n)/6+1):(3*n)/6}];
b=[varargin{(3*n)/6+1:(4*n)/6}];
E=[varargin{(4*n)/6+1:(5*n)/6}];
fl={varargin{((5*n)/6+1):n}};


%tensiones normales
Sn=@(F,A) F/A;
%tensiones flexoras
Sf=@(y,M,I,E) (y*M)/(E*I);
%tensiones de corte
Tc=@(b,V,I) (V*y^2)/(2*I);

%lo pongo todo en un mismo vector
Nn=n/6;
St=zeros(6,Nn);
for i=1:n/6

St(:,i)=[Sn(fl{i}(1),A(i));Tc(b(i),fl{i}(2),I(i));Sf(y(i),fl{i}(3),I(i),E(i));
        Sn(fl{i}(4),A(i));Tc(b(i),fl{i}(5),I(i));Sf(y(i),fl{i}(6),I(i),E(i))];
end

end

