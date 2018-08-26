% G2DFEAS
% General two dimensional Finite Element Analysis Script
% Requiere variables cargados
% Usar unidades CONSISTENTES: [N,mm,Mpa,mm^2,mm^4] [N,m,Pa,m^2,m^4] etc.
%   nod=[x1 y1;x2 y2... xn yn]
%   elenod=[N1 N2;N2 N3 ] conexiones entre nodos (elementos). Importante
%   tener bien definido el orden numerado de los elementos y los nodos!
%   Despues los resultados se basan en el orden que se tiene en elenod.
%   eletype=[1 2 2 1 3 4 1 2... 2 1 1];
%       1=barra
%       2=viga
%       3=viga rotulada en su nodo de inicio
%       4=viga rotulada en su nodo final o destino
%   apoyos_simples=[N1 N5 ... N8] Nodos apoyados en x e y
%   empotramientos=[] Nodos empotrados
% 
%       MATERIALES:
%   Ee=[200e9 210e9 ... 200e9] Modulos de young para cada elemento, se sigue el
%   orden de elenod.
%   Similarmente, se tiene que declarar Ae (area), Ie (2do momento de
%   inercia, he (altura/radio de viga, puede ser arbitrario para
%   barras), be (ancho/diametro para viga), Sye (tension de fluencia)
%   safetyfactor=4; factor seguridad deseado
% 
%       OPCIONES AVANZADAS (opcional):
%   CB=~~[1 1 0 1 1 ... 0  0  1]   Condiciones de borde (TRUE/1 anula grado de libertad)
%   Todos los nodos tienen 3 grados de liberad. Tomar esto en cuenta al
%   aplicar condiciones de borde. Las rotulas generan un grado de libertad
%   extra al final de CB, siguiendo el orden de elenod.
%   
%   graficar=true;  grafica estructura modelada
%   vigasinteresantes=[4:9] grafica esfuerzos sobre elementos 4 a 9 si
%   graficar=true.
% 
%       RESULTADOS:   
%   Muestra inmediatamente los resultados de tensiones sobre las barras.
%   dangerzone:
%   primera fila numera los elementos, segunda fila muestra el calculo:
%   (tension/tensionadm)*factorseguridad. Cuanto mas alto, mas
%   comprometido. En principio se desea que no exceda la unidad.
%   Segunda fila muestra peligro de pandeo para barras articuladas en los extremos.
%   
%       OUTPUT :
%   forzas: Devuelve estructura (fuerzas sobre cada elemento)
%   loskrotados: estructura con matriz rigidez locales.
%   
%   Patricio Whittingslow CC-BY-NC  2018 

%Actual code

verify(nod,elenod);
userelenod=elenod;
[nodeDofs,elenod,nudos,Le,phide,Ndof]=gennodedofs(nod,elenod,eletype);
[Nt,ndof]=size(nodeDofs);
% Ndof=ndof*Nt;
Ne=length(eletype);
T0b=[1 0;0 0;0 1;0 0];

sigmas=zeros(Ne,1);
taus=zeros(Ne,1);

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
elementos=zeros(Ne,6);
for i=1:Ne %ASSEMBLY
    switch eletype(i)
        case {1 , 11}
            klocal=Kb(Ee(i),Ae(i),Le(i));
            if eletype(i)==11
                klocal=klocal/2;
            end
            T=Tbu(phide(i));
            kbarrarotada=T*klocal*T';
            index=[nodeDofs(elenod(i,1),[1 2]) nodeDofs(elenod(i,2),[1 2])];
            kG(index,index)=kG(index,index)+kbarrarotada;
            klocalrotada=kbarrarotada;
            index=[nodeDofs(elenod(i,1),[1 2]) nodeDofs(elenod(i,2),[1 2])];
            elementos(i,[1 2 3 4 5 6])=...
                [index(1:2) nodeDofs(elenod(i,1),3) index(3:4) nodeDofs(elenod(i,2),3)];
            if nudos(elenod(i,1))==0
                CB(index(2)+1)=true;
            end
            if nudos(elenod(i,2))==0
                CB(index(4)+1)=true;
            end
        otherwise
            klocal=Kv(Ee(i),Ae(i),Ie(i),Le(i));
            if eletype(i)>11
                klocal=klocal/2;
            end
            T=Tvu(phide(i));
            klocalrotada=T'*klocal*T;
            index=[nodeDofs(elenod(i,1),:) nodeDofs(elenod(i,2),:)];
            elementos(i,:)=index;
            kG(index,index)=kG(index,index)+klocalrotada;
            
    end
     loskrotados=[loskrotados klocalrotada];%Guardo cada krotado
     losklocales=[losklocales klocal];
end

try
for i=1:length(apoyos_simples)
    n=apoyos_simples(i);
    CB([n*ndof-2 n*ndof-1])=true;
end
catch
end

try
for i=1:length(empotramientos)
    n=empotramientos(i);
    CB([n*ndof-2 n*ndof-1 n*ndof])=true;
end
catch
end

Kr=kG(~CB,~CB);
Rfulldof=zeros(Ndof,1);
Rfulldof(1:length(R))=R;
F=Rfulldof(~CB);

% F=R(~CB);
U=Kr\F; %OBTUVE DESPLAZANIETOS
D=zeros(Ndof,1);
D(~CB)=U;

forzasold={};
forzas={}; %Genero estructura con fuerzas sobre cada elemento
dangerzone=zeros(2,Ne);
for i = 1:Ne
    klocalr=loskrotados{i};
    klocal=losklocales{i};
    if eletype(i)==1
        ulocal=D(elementos(i,[1 2 4 5]));
        T=Tbu(phide(i));
        flocal=klocal*T'*ulocal;
        forzas=[forzas flocal];
    else
        ulocal=D(elementos(i,:));
        T=Tvu(phide(i));
        flocal=klocal*T*ulocal;%ROTACION DE FUERZAS A MI SISTEMA ELEMENTO
        forzas=[forzas flocal];
    end

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
        graficapoco(nod,userelenod,eletype,Ie);
    end
    if ~isempty(vigasinteresantes) && graficar      
        grafisuficiente2(Le(vigasinteresantes),forzas{vigasinteresantes})
    end
catch
    graficapoco(nod,userelenod,eletype,Ie);
end

% i=6;
% Calcutodo(Ae(i),Ie(i),ce(i),be(i),Ee(i),forzas{i});
try
    if showresults
        disp([1:Ne;dangerzone])
    end
catch
    disp([1:Ne;dangerzone])
end