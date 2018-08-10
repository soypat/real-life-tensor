
% TIEMPOS
timebs=0.001191; % Tiempo promedio para efectuar K\R con 10k nodos

clear all
pause(1);
E=30e6;
L=60;%in
A=2; %in^2
q=@(x) -10*x; %lb/in
N=10;
Ne=N-1;
h=L/Ne;
k=E/h*A;
tic
nodos=[0:h:L];
elementos=zeros(Ne,2);
for i=1:Ne
    elementos(i,1)=i;
    elementos(i,2)=i+1;
end

kG=sparse(N,N);

for i=1:Ne
    
    ke=k*[1 -1;-1 1];
    kG(elementos(i,:),elementos(i,:))=kG(elementos(i,:),elementos(i,:))+ke;
end
% CB=zeros(1,N);
R=zeros(N,1);
%Metodo parado sobre elementos
for i=1:Ne
    x=(i-1)*h;
    Atri=h/2*(q(x+h)-q(x));
    Acuadr=h*q(x);
    R(i)=R(i)+Acuadr/2+Atri/3;
    R(i+1)=Acuadr/2+Atri*2/3;
end

% for i=2:N-1 %Metodos conceptualmente incorrectos.
% %     Metodo Parado arriba de nodo. Promedio entre entrenodos
% %     x=h*i;
% %     R(i)=h/2*(q(x-h/2)+q(x+h/2));
%     
% %     Metodo entre-nodos;
% %     x=h*(i-1)-h/2;
% %     R(i)=h/2*(q(x)+q(x+h));
% %     son la misma cosa.
% end

P=q(L)*L/2;
% R(1)=q(h/2)*h/4;
R(end)=P;
Rr=R;
R=R(1:end-1);
K=kG(1:end-1,1:end-1);

% opts.POSDEF=true;
% U2= linsolve(K,R,opts);
U=K\R;
% upper=chol(K);
% Uc=(upper\(upper\eye(N-1))')*R;
% U=upper*(upper)
toc
% plot(nodos(1:Ne),U);
% 
% 
% xlabel('Posicion [in]');
% ylabel('Desplazamiento [in]')

%TENSIONES
D=zeros(Ne);
S=zeros(Ne);
for i=1:Ne-1
    D(i)=(U(i+1)-U(i))/h;
    S(i)=E*D(i);
end
nod=[h:h:L];
bar(nod,S,h);