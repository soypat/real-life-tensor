
% TIEMPOS
timebs=0.001191; % Tiempo promedio para efectuar K\R con 10k nodos

clear all
pause(1);
E=30e6;
L=60;%in
A=2; %in^2
q=@(x) -10*x; %lb/in
N=30000;
Ne=N-1;
h=L/Ne;
nodos=[0:h:L];
elementos=zeros(Ne,2);
for i=1:Ne
    elementos(i,1)=i;
    elementos(i,2)=i+1;
end

kG=sparse(N,N);

for i=1:Ne
    ke=E/h*A;
    ke=ke*[1 -1;-1 1];
    kG(elementos(i,:),elementos(i,:))=kG(elementos(i,:),elementos(i,:))+ke;
end
% CB=zeros(1,N);
R=z(N,1);
for i=2:N-1
    x=h*(i-1)-h/2;
    R(i)=h/2*(q(x)+q(x+h));
end

P=q(L)*L/2;
R(1)=q(h/2)*h/2;
R(end)=P;
R=R(2:end);
K=kG(1:end-1,1:end-1);

% opts.POSDEF=true;
% U2= linsolve(K,R,opts);
U=K\R;
% upper=chol(K);
% Uc=(upper\(upper\eye(N-1))')*R;
% U=upper*(upper)
toc
plot(nodos(1:Ne),U);

xlabel('Posicion [in]');
ylabel('Desplazamiento [in]')

    