
E=200e9;
b=.05;
h=.15;
A=b*h;
% problema 1

    

% Rigidez local elemento de viga
%Modele la cosa triangular como 3 vigas, dos de las cuales son rotuladas


    
    %BARRA
    

nod=[0 0;0 .2;.2 .2+.15;.4 .2+.15+.15;.2+.2+.1 .2+.15+.15+.15/2;.2+.2+.1, .2+.15+.15/2;
    .2+.15 .2;.4 -.1;.2 .2];
N=size(nod,1);
% scatter(nod(:,1),nod(:,2))
enod=[1 2;2 3;3 4;4 5;4 6;6 5;7 3;7 9;9 3;8 9];
eletype=[1 2 2 2 2 1, 3 2 3, 4];% 1 barra, 2 viga, 3 rotula final, 4 piston
Ne=size(enod,1);

%assemble bars
fnod=@(n) [n*3-2 n*3-1 n*3];
fnodb=@(n) [n*3-2 n*3-1];
% dex1=[1];
Le=zeros(Ne,1);
phide=zeros(Ne,1);
Ndof=N*3+1;%3 grados libertad para todos, despues los desaparezco con CB
IN=zeros(Ne,6);%es la matriz elementos
kG=zeros(Ndof);
Ts={};
for i=1:Ne
    ns=enod(i,1);
    ne=enod(i,2);
    lx=nod(ne,1)-nod(ns,1);
    ly=nod(ne,2)-nod(ns,2);
    Le(i)=sqrt(lx^2+ly^2);
    phide(i)=atan2d(ly,lx);
    c=cosd(phide(i));
    s=sind(phide(i));
    T=[c s 0 0 0 0;-s c 0 0 0 0;0 0 1 0 0 0;0 0 0 c s 0;0 0 0 -s c 0;0 0 0 0 0 1];
    Tb=[c 0;s 0;0 c;0 s];
    Ts=[Ts T];
    klocs={};
    switch eletype(i)
        case 2
            Iz=b*h^3/12; %VIGA
            nu=.3;
            X = E*A/Le(i);
            Y1 = 12*E*Iz/Le(i)^3;
            Y2 = 6*E*Iz/Le(i)^2;
            Y3 = 4*E*Iz/Le(i)^1;
            Y4 = 2*E*Iz/Le(i)^1;
            Kdiag = diag([X Y1 Y3 X Y1 Y3]);
            Kp = [0     0       0       -X      0       0
                  0     0       Y2      0       -Y1     Y2
                  0     0       0       0       -Y2     Y4
                  0     0       0       0       0       0
                  0     0       0       0       0       -Y2
                  0     0       0       0       0       0];
            Kel = Kp + Kp' + Kdiag;
            
            if i==2
                klocalBC=Kel;
            end
            if i==3
                klocalCD=Kel;
            end
            klocs=[klocs Kel];
            index=[fnod(ns) fnod(ne)];
            krot=T'*Kel*T;
            IN(i,:)=index;
 %%   INCISO 2 (con.)
            if i==8
                krot=10*krot;%PARTE DEL INCISO 2. 
            end
        case 1 %Calculo de barras.
      %%      %%INCISO 1 
            %Le agrego un inciso al codigo para tomar mi elemento 6 con mas
            %rigidez pues no forma parte de la estructura que quiero
            %analizar. 10 veces mas rigidez para ser exactos.
            if i==6
                Kel = E*A*10/Le(i)*[1 -1;-1 1];
                index=[fnodb(ns) fnodb(ne)];
                IN(i,[1 2 4 5])=index;
                krot=Tb*Kel*Tb';
                kG(index,index)=kG(index,index)+krot;
                continue
          
            end%%Hasta aqui agregue

            Kel = E*A/Le(i)*[1 -1;-1 1];
            index=[fnodb(ns) fnodb(ne)];
            IN(i,[1 2 4 5])=index;
            krot=Tb*Kel*Tb';
            
            
        case 4%Calculo de piston
            Ap=.05^2*pi/4;
            Kel = E*Ap/Le(i)*[1 -1;-1 1];
            kp=Kel;
            index=[fnodb(ns) fnodb(ne)];
            krot=Tb*Kel*Tb';
            Tp=Tb;
            IN(i,[1 2 4 5])=index;
        case 3
%%            %INCISO 2
            %Tambien creo que deberia agregarle mas rigidez a los elementos
            %que forman parte del eslabon triangular, pues no son
            %conformados por secciones como los de las vigas. todas las
            %rotulas son del eslabon asi que generalizo
            Iz=b*h^3/12; %Rotulas
            nu=.3;
            X = 5*E*A/Le(i);
            Y1 = 12*E*Iz/Le(i)^3;
            Y2 = 6*E*Iz/Le(i)^2;
            Y3 = 4*E*Iz/Le(i)^1;
            Y4 = 2*E*Iz/Le(i)^1;
            Kdiag = diag([X Y1 Y3 X Y1 Y3]);
            %PARTE DEL INCISO 2 -> Multiplico las vigas por 10
            Kp = [0     0       0       -X      0       0
                  0     0       Y2      0       -Y1     Y2
                  0     0       0       0       -Y2     Y4
                  0     0       0       0       0       0
                  0     0       0       0       0       -Y2
                  0     0       0       0       0       0];
            Kel = Kp + Kp' + Kdiag;
            Kel=Kel*10; %ACA!
            
            krot=T'*Kel*T;
            index=[fnod(ns) fnod(ne)];
            index(6)=N*3+1;%mando la rigidez a la rotula
            IN(i,:)=index;
    end
        kG(index,index)=kG(index,index)+krot; %Aqui acoplo mi local
end
    %% INCISO 3 -> me olvide de no incluir el giro en nodo 8.
disp(kG([2:23 25:end],[2:23 25:end])) %Mi matriz global. Ojo, todos los nodos les di 3 grados de libertad. se puede indexear los k rigidez para que te quede mas chico, en particular para u1 que da cero
    %%
CB=false(Ndof,1);
CB([fnodb(8) fnodb(7)])=true; %Apoyos simples condicionados
CB([fnod(1) fnod(8)])=true;%Apoyos simples condicionados con barra
CB
R=zeros(Ndof,1);
P1=50;  %N y Nm
P2=100;
M1=8.333333;
M2=12.5;

R(fnod(5))=[0 -(P1+P2)/2 (M1+M2)/2];%divido las fuerzas equivalently 
R(fnod(6))=[0 -(P1+P2)/2 (M1+M2)/2];

F=R(~CB);
Kr=kG(~CB,~CB);

D=zeros(Ndof,1);

U=inv(Kr)*F;

D(~CB)=U;

dezp_plat=D([fnod(5) fnod(6)]);
giros=dezp_plat(3);% GIROS 1b

ulocalB=D([fnod(2) fnod(3)])
flocal=klocalBC*Ts{2}*ulocalB
Mbc=(flocal(3)+flocal(6))*.5;%Se que es lineal, la aproximo como tal
N=-flocal(1);
V=flocal(2); %EN PASCALES
sigmaBC=N/A+abs(Mbc)*h/2/Iz

tauBC=3/2*V/A %secciones cuadradas

up=D([fnodb(8) fnodb(9)]);

flocalp=kp*Tp'*up ;%No estaria llegando a calcular estas
N=-flocalp(1);
sigmaPiston=N/Ap

%% Anexo.
%me falto calcuar las fuerzas sobre el Bulon del eslabon
%triangular.----------------------------------------------------------
flocalele3=klocalCD*Ts{3}*D([fnod(3) fnod(4)]);
fuerzas_xprimayprima=flocalele3([4 5])+flocal([1 2]); %fuerzas primadas
FuerzaTotalEnElBulonSiTeInteresa=norm(fuerzas_xprimayprima)
%Y asi creo que terminaba el ultimo item que me faltaba.

