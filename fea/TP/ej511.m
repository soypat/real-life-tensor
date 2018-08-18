A=8;
I=800;
E=30e6;
Ae=A*ones(Ne,1);
Ee=E*ones(Ne,1);
Ie=I*ones(Ne,1);
ndof=3;
nodes=[0 0;
    40 0;
    70 40;
    20 40];
nodes=nodes*12;
elenod=[1 4;
        2 4;
        3 4];
eletype=[1 1 1]
%INPUTS END
[Ne,~]=size(elenod);
[N,~]=size(nodes);
elementos=genelementos(elenod);
Le=zeros(Ne,1);
phide=zeros(Ne,1);
Ndof=N*ndof;
lx=zeros(Ne,1);
ly=zeros(Ne,1);

for i = 1:Ne
    nodestart=elenod(i,1);
    nodeend=elenod(i,2);
    lx(i)=nodes(nodeend,1)-nodes(nodestart,1);
    ly(i)=nodes(nodeend,2)-nodes(nodestart,2);
    Le(i)=sqrt(lx(i)^2+ly(i)^2);
    phide(i)=atan2d(ly(i),lx(i));%angulo en degrees
end

kG=zeros(Ndof);
losklocales={};
loskrotados={};

for i=1:Ne
    klocal=Kv(Ee(i),Ae(i),Ie(i),Le(i));
    losklocales=[losklocales klocal];
    T=Tvu(phide(i));
    klocalrotada=T'*klocal*T;
    loskrotados=[loskrotados klocalrotada];
    kG(elementos(i,:),elementos(i,:))=kG(elementos(i,:),elementos(i,:))+klocalrotada;
end



CB=false(Ndof,1);
CB([1 2 3 4 5 6 7 8 9])=true;

K=kG(~CB,~CB); %Generada Bien!



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
