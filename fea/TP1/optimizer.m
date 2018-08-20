Ao=Ae;
co=ce;
bo=be;
ho=he;
Io=Ie;
No=4;
kG=zeros(Ndof);
optielementos=[];%[1 3 4 5 9];%[6 7 8];  %[1 3 4 5 6 7 8 9];
safetyfactor=6;
affinity=1.2; %que tan fino queremos la aproximacion al optimo. Mas grande, mas fino.
No=1000;
for j=1:1
    sigmas=zeros(Ne,1);
    taus=zeros(Ne,1);
    for i=1:Ne %ASSEMBLY
        switch eletype(i)
            case 1
                klocalrotada=Kb(Ee(i),Ao(i),Le(i),phide(i));
                CB(elementos(i,ndof))=true;
                CB(elementos(i,2*ndof))=true;
            case 2
                klocal=Kv(Ee(i),Ao(i),Io(i),Le(i));
                T=Tvu(phide(i));
                klocalrotada=T'*klocal*T;
            case 3
                klocal=vigabisagrada(Ee(i),Ao(i),Io(i),Le(i),1);
                T=Tvu(phide(i));
                klocalrotada=T'*klocal*T;
            case 4
                klocal=vigabisagrada(Ee(i),Ao(i),Io(i),Le(i),0);
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
        if eletype(i)>1
            [sig, tau]=getvigatensions(bo(i),ho(i),phide(i),flocal);
        else
            sig=getbartensions(bo(i),phide(i),flocal);
            tau=0;
        end
        sigmas(i)=sig;
        taus(i)=tau;
        dangerzone(i)=abs(sig)/Sye(i)*safetyfactor;
    end
    %COMIENZA OPTIMIZACION
    for i=optielementos
        dz=dangerzone(i);
%         if dz<0
%             continue
%         end
        slopedy=pi/(24*affinity);
        if abs(dangerzone(i))<1 %&& dangerzone(i)>0 %Entramos a rutina de agregarle material
            roughnessindex=(1-abs(dangerzone(i)))*affinity; %Nos dice que tan lejos estamos de nuestra meta
            ho(i)=ho(i)-min(roughnessindex,slopedy);
            if bo(i)>1
                bo(i)=bo(i)-min(roughnessindex,slopedy);
            end
        elseif abs(dangerzone(i))>1 %entramos a rutina de Sacarle material
            roughnessindex=(dangerzone(i)-1)*affinity; %Nos dice que tan lejos estamos de nuestra meta
            ho(i)=ho(i)+min(roughnessindex,slopedy);
            if bo(i)>1
                bo(i)=bo(i)+min(roughnessindex,slopedy);
            end
            
        end
        %Update Areas y Momentos de inercia
        if eletype(i)>1 %Para vigas
            Ao(i)=bo(i)*ho(i);
            Io(i)=bo(i)*ho(i)^3/12;
        else
            Ao(i)=bo(i)^2*pi/4;%Las barras son redondas
            Io(i)=bo(i)^4*pi/64;
        end
    end
end