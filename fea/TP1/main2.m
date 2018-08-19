%Main v2
%G2D-FEAS
%General two dimensional Finite Element Analysis Script
%Preparation:
enunciado

%Actual code
Le=zeros(Ne,1);
phide=zeros(Ne,1);
Ndof=N*ndof;

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
% GrafitodoF(Le(6),Le(7),Le(8),forzas{6},forzas{7},forzas{8});
graficapoco(nod,elenod,eletype,Ie);

grafisuficiente(Le(6:9),forzas{6:9})

i=6;
Sy=Calcutodo(Ae(i),Ie(i),ce(i),be(i),Ee(i),forzas{i});