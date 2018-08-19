kG=zeros(Ndof);
for i=1:Ne
    I
    
end
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

try
    if ~isempty(vigasinteresantes)
        grafisuficiente(Le(vigasinteresantes),phide(vigasinteresantes),forzas{vigasinteresantes})
    end
catch
end