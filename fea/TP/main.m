%Retroexcavadora
%barras
%
%Tb

%T(phi)
%viga(E,I,L)

%vigota(E,A,I,L)
%Tvu(phi)
hingotaL=@(E,A, I, L) (3*E*I/L^3)*[A*L^2/(3*I) 0 0 -A*L^2/(3*I) 0 0;
    0 1 L 0 -1 0;
    0 L L^2 0 -L 0;
    -A*L^2/(3*I) 0 0 A*L^2/(3*I) 0 0;
    0 -1 -L 0 1 0;
    0 0 0 0 0 0];
hingotaR=@(E,A,I,L) (3*E*I/L^3)*[A*L^2/(3*I) 0 0 -A*L^2/(3*I) 0 0;
    0 1 0 0 -1 L;
    0 0 0 0 0 0;
    -A*L^2/(3*I) 0 0 A*L^2/(3*I) 0 0;
    0 -1 0 0 1 -L;
    0 L 0 0 -L L^2];
P=3200;
E=30e6;
h_nom=4;%in
b_nom=3;
Fa=-P*15/16;
Fb=-Fa;
elementos=[1 2 3 4 5 6;
    1 2 3 10 11 12;
    7 8 9 10 11 12;
    4 5 6 7 8 9;
    7 8 9 13 14 15;
    10 11 12 13 14 15;
    13 14 15 16 17 18;
    16 17 18 19 20 21;
    16 17 18 22 23 24];

L1=20;
L2=16;
L3 = dist(20,8);
L4=8;
L5=dist(15,5);
L6=dist(7,15);
L7=dist(16,12);
L8=dist(6,24);
L9=dist(30,30);
phi1=0;
phi2=-90;
phi3=180+atand(8/20);
phi4=-90;
phi5=180+atand(15/5);
phi6=-atand(7/15);
phi7=-atand(12/16);
phi8=-atand(21/6);
phi9=-atand(1);
Iz=@(b,h) b*h^3/12;
Ar=@(b,h) b*h;
barras=[1 2 3 9];

h=h_nom;
b=b_nom;
Ne=9;
ndof=3;
N=8;
Ndof=ndof*N;
kG=zeros(Ndof);
I=Iz(h,h);
losKlocos={};
quieroanalizarlasmatrices=[6 7 8 9]; %QUIERO ANALIZAR LAS MATRICES = qalm
qalm=quieroanalizarlasmatrices; %Quiero Analizar Las Matrices   qalm
for i=1:9
    L=eval(sprintf('L%0.0f',i));
    A=Ar(h,h);
    phi=eval(sprintf('phi%0.0f',i));
    if sum(ismember(barras,i))>0%si es barra corre
         klocalrotada=Kb(E,A,L,phi);
         kG(elementos(i,:),elementos(i,:))=kG(elementos(i,:),elementos(i,:))+klocalrotada;
    elseif i==2
        Eeffectiva=E*10;
        klocalrotada=Kb(E,A,L,phi);
        kG(elementos(i,:),elementos(i,:))=kG(elementos(i,:),elementos(i,:))+klocalrotada;
    elseif i==5
        klocal=hingota(E,A,I,L);
        T=Tvu(phid);
        klocalrotada=T'*klocal*T;
        kG(elementos(i,:),elementos(i,:))=kG(elementos(i,:),elementos(i,:))+klocalrotada;
    else
        I=Iz(h,h);
        klocal=Kv(E,A,I,L);
        T=Tvu(phi);
        klocalrotada=T'*klocal*T;
        kG(elementos(i,:),elementos(i,:))=kG(elementos(i,:),elementos(i,:))+klocalrotada;

    end
    if sum(ismember(qalm,i))>0%guardo los valores klocales para vigas que Quiero Analizar
        losKlocos=[losKlocos klocal];
    end
end


CB=false(Ndof,1);

% CB(elementos(6,3))=true; %Viga

CB(elementos(8,4:5))=true;
CB(elementos(9,4:6))=true;
CB(elementos(1,6))=true;
CB(elementos(1,3))=true;

CB(elementos(2,6))=true;
CB(elementos(2,3))=true;
CB(elementos(3,6))=true;
CB(elementos(3,3))=true;

R=zeros(Ndof,1);
R(elementos(1,1:2))=[Fa -P/2];
R(elementos(2,4:5))=[Fb -P/2];

F=R(~CB);
K=kG(~CB,~CB);

U=K\F;
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
forzas={};
for i = 1:length(losKlocos)
    klocal=losKlocos{i};
    ulocal=D(elementos(qalm(i),:));
    
    flocal=klocal*ulocal;
    forzas=[forzas flocal];
end




% 
%     elseif i==6
%         kcita=hingotaR(E,A,I,L);
%         T=Tvu(phid);
%         klocal=T'*kcita*T;
%         kG(elementos(i,:),elementos(i,:))=kG(elementos(i,:),elementos(i,:))+klocal;