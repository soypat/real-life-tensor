%Main v2
%G2D-FEAS
%General two dimensional Finite Element Analysis Script

%Preparation:
% enunciado
enunciado
% vigasinteresantes=[6:8];
% test
ndof=3;
empotramientos=[];
%Actual code
Le=zeros(Ne,1);
phide=zeros(Ne,1);
hinges=length(eletype(eletype>2));
Ndof=N*ndof+hinges;
T0b=[1 0;0 0;0 1;0 0];
verify(nod,elenod);
sigmas=zeros(Ne,1);
taus=zeros(Ne,1);
nudos=zeros(Ndof,1);
for i = 1:Ne
    nodestart=elenod(i,1);
    nodeend=elenod(i,2);
    lx=nod(nodeend,1)-nod(nodestart,1);
    ly=nod(nodeend,2)-nod(nodestart,2);
    Le(i)=sqrt(lx^2+ly^2);
    phide(i)=atan2d(ly,lx);%angulo en degrees
    switch eletype(i)
        case 2
            nudos(nodestart*ndof)=nudos(nodestart*ndof)+1;
            nudos(nodeend*ndof)=nudos(nodeend*ndof)+1;
        case 3
            nudos(nodeend*ndof)=nudos(nodeend*ndof)+1;
        case 4
            nudos(nodestart*ndof)=nudos(nodestart*ndof)+1;
    end
end

kG=zeros(Ndof);
losklocales={};
loskrotados={};

try %Generacón del vector CB. Asegura correcta generación si el usuario tira fruta
    [a, b]=size(CB);
    if (a==1 && b==Ndof)
        CB=CB';
    elseif (a~=Ndof || ~islogical(CB))
        warning('Condiciones de borde no tienen tamaño requerido o no es matriz logica. Usar false().\nImponiendo condiciones de borde sugeridas...')
        CB=false(Ndof,1);
    end
catch e
    warning('Imponiendo condiciones de borde sugeridas debido a error relacionado con condiciones de borde...\nMensaje:%s',e.message)
    CB=false(Ndof,1);
end

for i=1:Ne %ASSEMBLY
    switch eletype(i)
        case 1
            kbarra=Kb(Ee(i),Ae(i),Le(i));
            klocal=zeros(6);
            kbarra6=T0b*kbarra*T0b';
            klocal([1 2 4 5],[1 2 4 5])=kbarra6;
            T=Tbu(phide(i));
            kbarrarotada=T*kbarra*T';
            klocalrotada=zeros(6);
            klocalrotada([1 2 4 5],[1 2 4 5])=kbarrarotada;
            if nudos(elementos(i,ndof))==0
                CB(elementos(i,ndof))=true;
            end
            if nudos(elementos(i,2*ndof))==0
                CB(elementos(i,2*ndof))=true;
            end
            kG(elementos(i,:),elementos(i,:))=kG(elementos(i,:),elementos(i,:))+klocalrotada;
        case 2
            klocal=Kv(Ee(i),Ae(i),Ie(i),Le(i));
            T=Tvu(phide(i));
            klocalrotada=T'*klocal*T;
            kG(elementos(i,:),elementos(i,:))=kG(elementos(i,:),elementos(i,:))+klocalrotada;
        case 3
            klocal=vigabisagrada(Ee(i),Ae(i),Ie(i),Le(i),1);
            T=Tvu(phide(i));
            klocalrotada=T'*klocal*T;
            kG(elementos(i,:),elementos(i,:))=kG(elementos(i,:),elementos(i,:))+klocalrotada;
                        
        case 4
%             klocal=vigabisagrada(Ee(i),Ae(i),Ie(i),Le(i),0);
            klocal=Kv(Ee(i),Ae(i),Ie(i),Le(i));
            Tvu(phide(i));
            klocalrotada=T'*klocal*T;
            kG(elementos(i,:),elementos(i,:))=kG(elementos(i,:),elementos(i,:))+klocalrotada;
    end
     kG(elementos(i,:),elementos(i,:))=kG(elementos(i,:),elementos(i,:))+klocalrotada;
     loskrotados=[loskrotados klocalrotada];%Guardo cada krotado
     losklocales=[losklocales klocal];
end



for i=1:length(apoyos_simples)
    n=apoyos_simples(i);
    CB([n*ndof-2 n*ndof-1])=true;
end

CB(6)=false;
Kr=kG(~CB,~CB);
F=R(~CB);
U=Kr\F; %OBTUVE DESPLAZANIETOS

D=regendesplazamientos(U,CB);


forzasold={};
forzas={}; %Genero estructura con fuerzas sobre cada elemento
dangerzone=zeros(2,Ne);
for i = 1:Ne
    klocalr=loskrotados{i};
    klocal=losklocales{i};
    ulocal=D(elementos(i,:));
    T=Tvu(phide(i));
    flocalold=klocalr*ulocal;
    flocal=klocal*T*ulocal;%ROTACION DE FUERZAS A MI SISTEMA ELEMENTO
    forzasold=[forzasold flocalold];
    forzas=[forzas flocal];
    if eletype(i)>1
        [sig, tau, N, M]=getvigatensions(be(i),he(i),flocal);
    else
        [sig, N]=getbartensions(be(i),flocal);
        tau=0;
        M=0;
    end
    if N<0
        Pcrit=pi^2*Ee(i)*Ie(i)/Le(i)^2;
    else
        Pcrit=NaN;
    end
    sigmas(i)=sig;
    taus(i)=tau;
    dangerzone([1,2],i)=[sig/Sye(i)*safetyfactor;-N/Pcrit*safetyfactor];
end
% GrafitodoF(Le(6),Le(7),Le(8),forzas{6},forzas{7},forzas{8});


try
    if graficar
        graficapoco(nod,elenod,eletype,Ie);
    end
    if ~isempty(vigasinteresantes) && graficar      
        grafisuficiente2(Le(vigasinteresantes),forzas{vigasinteresantes})
    end
catch
end

% i=6;
% Calcutodo(Ae(i),Ie(i),ce(i),be(i),Ee(i),forzas{i});
disp([1:Ne;dangerzone])