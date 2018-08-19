%Main v2
nod=[0 0;20 0; 20 -8;0 -16;15 -(15+8);%nodos 1 2 3 4 5
    (20+16) -(8+15+12);(16+20+24) -(6+12+15+8);(20+16+24+6) -(8+15+12+6+24)]; 
elenod=[1 2;1 4;3 4;2 3;3 5;4 5;5 6;6 7;6 8];
ndof=3;
[N,~]=size(nod);
Ndof=N*ndof;

[Ne,~]=size(elenod);
Ee=30e6*ones(Ne,1);
Ae=8*ones(Ne,1);
Ie=800*ones(Ne,1);
Ie(2)=Ie(2)*4;
Ie(6)=Ie(6)*8;
Ie(7)=Ie(7)*8;
P=3200;
Fa=-P*15/16;
Fb=-Fa;



eletype=[1 1 1 3 4 2 2 2 1];
%Types
%1=barra
%2=viga comun
%3=viga con bisagra al comienzo
%4=viga con bisagra al final


elementos=genelementos(elenod);
Le=zeros(Ne,1);
phide=zeros(Ne,1);
Ndof=N*ndof;
R=zeros(Ndof,1);
R(elementos(1,1:2))=[Fa -P/2];
R(elementos(2,4:5))=[Fb -P/2];
for i = 1:Ne
    nodestart=elenod(i,1);
    nodeend=elenod(i,2);
    lx=nod(nodeend,1)-nod(nodestart,1);
    ly=nod(nodeend,2)-nod(nodestart,2);
    Le(i)=sqrt(lx^2+ly^2);
    phide(i)=atan2d(ly,lx);%angulo en degrees
end

kG=zeros(Ndof);
losklocales={};
loskrotados={};
CB=false(Ndof,1);

for i=1:Ne %ASSEMBLY
    switch eletype(i)
        case 1
            klocalrotada=Kb(Ee(i),Ae(i),Le(i),phide(i));
            CB(elementos(i,ndof))=true;
            CB(elementos(i,2*ndof))=true;
        case 2
            klocal=Kv(Ee(i),Ae(i),Ie(i),Le(i));
            T=Tvu(phide(i));
            klocalrotada=T'*klocal*T;
        case 3
            klocal=vigabisagrada(Ee(i),Ae(i),Ie(i),Le(i),1);
            T=Tvu(phide(i));
            klocalrotada=T'*klocal*T;
        case 4
            klocal=vigabisagrada(Ee(i),Ae(i),Ie(i),Le(i),0);
            Tvu(phide(i));
            klocalrotada=T'*klocal*T;
    end
     kG(elementos(i,:),elementos(i,:))=kG(elementos(i,:),elementos(i,:))+klocalrotada;
     loskrotados=[loskrotados klocalrotada];%Guardo cada krotado
end

apoyos_simples=[7 8];

for i=1:length(apoyos_simples)
    n=apoyos_simples(i);
    CB([n*ndof-2 n*ndof-1])=true;
end
% clf
Kr=kG(~CB,~CB);
F=R(~CB);
U=Kr\F; %OBTUVE DESPLAZANIETOS

D=regendesplazamientos(U,CB);


forzas={}; %Genero estructura con fuerzas sobre cada elemento
for i = 1:Ne
    klocal=loskrotados{i};
    ulocal=D(elementos(i,:));
    flocal=klocal*ulocal;
    forzas=[forzas flocal];
end

graficapoco(nod,elenod,eletype,Ie);


% function [D]=regendesplazamientos(U,CB)
% Ndof=length(CB);
% CBneg=~CB;
% D=zeros(Ndof,1);
% Ucount=0;
% %REgenera el vector U con 0 en desplazamientos restringidos
% for i=1:length(CB)
%     if CBneg(i)
%         Ucount=Ucount+1;
%         D(i)=U(Ucount);
%     else
%         continue
%     end
% end
% end

% function [] = graficapoco(nod,elenod)
% clf
% [Ne,~]=size(elenod);
% scatter(nod(:,1),nod(:,2),'.k')
% for i = 1:Ne
%     hold on
%     xv=[nod(elenod(i,1),1) nod(elenod(i,2),1)];
%     yv=[nod(elenod(i,1),2) nod(elenod(i,2),2)];
%     line(xv,yv);
% end
% 
% end