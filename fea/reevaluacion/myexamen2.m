b=.1;
h=.3;
A1=(15^2-10^2)*pi/4;
A2=(10^2)*pi/4;
A3=5^2;
Iz1=(15^4-10^4)*pi/64;
Iz2=10^4*pi/64;
Iz3=5*5^3/12;

E=200e3; %Mpa
nod=[0 0;100 0;125 0;200 0; 75 -25;125 -25;75 -75;125 -75;100 -100];
elenod=[1 2;2 3;3 4;2 5;2 6;5 7;5 8;6 8;8 9;7 9];
Ne=size(elenod,1);
N=size(nod,1);
Le=zeros(Ne,1);
phide=zeros(Ne,1);
fnod=@(n) [n*3-2 n*3-1 n*3];
fnodb=@(n) [n*3-2 n*3-1];
Tcell=cell(Ne,1);
kG=zeros(N*3);
kloc=cell(Ne,1);
for e=1:Ne
    ns=elenod(e,1);
    ne=elenod(e,2);
    lx=nod(ne,1)-nod(ns,1);
    ly=nod(ne,2)-nod(ns,2);
    
    px=nod(ne,:)-nod(ns,:);
    Le(e)=norm(px);
    pn=px/norm(px);
    C=[pn(1) pn(2) 0;-pn(2) pn(1) 0;0 0 1];
    T=blkdiag(C,C);
    Tcell{e}=T;
    switch e
        case {1 2}
            Iz=(15^4-10^4)*pi/64;
            A=(15^2-10^2)*pi/4;
            
%             Iz=(15^4)*pi/64;%Caso carga inicial
%             A=(15^2)*pi/4;

        case 3
            A=(10^2)*pi/4;
            Iz=10^4*pi/64;

        otherwise
            A=5^2;
            klocal= E*A/Le(e)*[1 -1;-1 1];
            index=[fnodb(ns) fnodb(ne)];
            Tb=[pn(1) pn(2) 0 0;0 0 pn(1) pn(2)];
            
            krotada=Tb'*klocal*Tb;
            kG(index,index)=kG(index,index)+krotada;
            kloc{e}=klocal;
            continue
    end
    Kel=Kv(E,A,Iz,Le(e));
    krotada=T'*Kel*T;
    index=[fnod(ns) fnod(ne)];
    %%Debugging problema mal condicionamiento
    if e==3
        kG(index,index)=kG(index,index)+krotada*2;
    else
        kG(index,index)=kG(index,index)+krotada;
    end
end 

CB=false(3*N,1);

% CB(fnod(4))=[false true false]; %WARNING Esta mal condicionada mi matriz. para solucionarlo lo unico que encontre fue aplicar condicion
%de borde en mi nodo C (4) en direccion Y. Aun asi no se me arregla el
%problema (obviamente). Voy a tener que encontrar la raiz del problema en
%casa

CB(fnod(5))=[false false true]; %condiciona barras
CB(fnod(6))=[false false true];
CB(fnod(7))=[false false true];
CB(fnod(8))=[false false true];

CB(fnod(9))=[true true true];%Apoyo
%%TENGO UNA MATRIZ MAL CONDICIONADA (K reducida)
CB([1 2])=true; %apoyo A
kR=kG(~CB,~CB);

R=zeros(3*N,1);
R(fnod(3))=[0 -2000 0];
F=R(~CB);
U=kR\F;
%%Me da que la viga se desplaza a~nos luz. Un error mio. 
D=zeros(3*N,1);
D(~CB)=U;
D(fnod(3))%
% %% NExt part
R2=zeros(N*3,1);
R2(fnod(4))=[1000 -1000 0]*sind(45);
F2=R2(~CB);
U2=kR\F2;
D2=zeros(3*N,1);
D2(~CB)=U2;
nodo4=fnod(4);
%% desplazamientos y giros 2 a
D2(fnod(4))

%% viga AB piston
Dab=D2([fnod(1) fnod(2)]);
Iz=(15^4-10^4)*pi/64;
A=(15^2-10^2)*pi/4;
L=Le(1);
x=L/2;
Mz=E*Iz*( (-6/L^2 + 12*x/L^3)*Dab(2)+(-4/L+6*x/L^2)*Dab(3)+(6/L^2-12*x/L^3)*Dab(5)+(-2/L+6*x/L^2)*Dab(6) )
Nx=Dab(1)*-1/L+Dab(4)*1/L
sigma=Nx/A+(Mz*(15/2)/Iz)
