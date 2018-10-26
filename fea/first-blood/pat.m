%In class

E=210e3;
P=1000;%N
L=4000;
N=5;
Ne=N-1;
h=L/Ne;%Paso
% A=@(x) 25*(1+x*3/4000);
A=@(x) 25*(L-x)/L+100*x/L;
nodos=[0:h:L];
elementos=zeros(Ne,2);
for i=1:Ne
    elementos(i,1)=i;
    elementos(i,2)=i+1;
end

kG=sparse(N,N);

for i=1:Ne
    ke=E/h*(A(h*i));
    ke=ke*[1 -1;-1 1];
    kG(elementos(i,:),elementos(i,:))=kG(elementos(i,:),elementos(i,:))+ke;
end
CB=ones(1,N);
CB(end)=0;
% K=kG()
R=zeros(Ne,1);
R(1)=-P;

CB2=logical(CB);

K=kG(CB2,CB2);
tic;
% K=kG(2:end,2:end);
U=K\R;

plot(nodos(1:Ne),U);
toc
clear U
pause(1);

plot(nodos(1:Ne),U);


%First try
% A1=37.5;
% A2=62.5;
% A3=87.5;
% A1=50;
% A2=75;
% A3=100;
% E=210e6; %MPa
% Le=4/3;%m
% k1=E*A1/Le;
% k2=E*A2/Le;
% k3=E*A3/Le;
% P=1000;
% 
% K=[k1 -k1 0 0;
%     -k1 k1+k2 -k2 0;
%     0 -k2 k2+k3 -k3;
%     0 0 -k3 k3];
% Kutil=K(2:end,2:end);
% Kinv=Kutil^-1;
% R=[0;0;P];
% U=Kinv*R;

%Generalizeations
% L=4000;%mm
% N=4000;
% Ne=N-1;
% Le=L/N;
% K=zeros(N);
% A=@(x) 25*(1+x*3/4);
% for i=1:N
%     for j=1:N
%         if abs(i-j)>1
%             K(i,j)=0;
%             continue;
%         end
%         if i==j
%             kij=E/Le*(A(Le*i-Le/2));
%             K(i,j)=kij;
%             continue;
%         else
%             l=min(i,j);
%             kij=E/Le*(A(Le*l-Le/2));
%             K(i,j)=-kij;
%         end
%     end
% end
% Kutil=K(2:end,2:end);
% R=zeros(N-1,1);
% R(end)=P;
% U=Kutil\R;
% plot(U)

