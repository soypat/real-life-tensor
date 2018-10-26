%ejemplo 43
nod=[4 4 3;%1
       0 4 0;%2
       0 4 6;%3
       4 0 3;%4
       8 -1 1];%5
enod=[1 2;
        1 3;
        1 4;
        1 5];

fnod=@(n) [n*6-5 n*6-4 n*6-3 n*6-2 n*6-1 n*6];
[Ne,~]=size(enod);
[N,~]=size(nod);
Le=zeros(Ne,1);
% C=zeros(Ne,1);


Ts={};
kG=zeros(6*N);
Ee=210e9*ones(Ne,1);%Pa
Ae=10e-4*ones(Ne,1);%m
for i=1:Ne
    ns=enod(i,1);
    ne=enod(i,2);
    lx=nod(ne,1)-nod(ns,1);
    ly=nod(ne,2)-nod(ns,2);
    lz=nod(ne,3)-nod(ns,3);
    Le(i)=sqrt(lx^2+ly^2+lz^2);%Long.
    px=[lx ly lz]/Le(i);
    vy=[0 1 5];
    pz=cross(px,vy)/norm(cross(px,vy));
    py=cross(pz,px);
    C=[px;py;pz;];%Para barras no importa orientacion azimuth (acimutal)
    T=blkdiag(C,C,C,C);
    Ts=[Ts T];
    klocal=Kvuw(Ee(i),Ae(i),0,0,0,0,Le(i));
    krotada=T'*klocal*T;
    kG([fnod(enod(i,1)) fnod(enod(i,2))],[fnod(enod(i,1)) fnod(enod(i,2))])=...
      kG([fnod(enod(i,1)) fnod(enod(i,2))],[fnod(enod(i,1)) fnod(enod(i,2))])+krotada;
end
CB=false(6*N,1);
CB([fnod(2) fnod(3) fnod(4) fnod(5)])=true;
CB(fnod(1))=~~[0 0 0 1 1 1];
Kr=kG(~CB,~CB);
R=zeros(6*N);
R(fnod(1))=[0 -10e3 0 0 0 0];
F=R(~CB);
U=Kr\F;
disp(U)
 %Nos aseguramos que no haya forma de obtener U accidentalmente 
%Ahora imponemos desplazamiento inicial
% v1=-1.517731e-04;
% Dc=zeros(N*6,1);
% Dc(fnod(1))=[0 v1 0 0 0 0];
% CB(fnod(1))=~~[0 1 0 1 1 1];
% R=zeros(N*6,1); %si hay fuerzas ademas de la desconocida las agregamos
% Dx=kG(~CB,~CB)\(R(~CB)-kG(~CB,CB)*Dc(CB));
% Rx=kG(CB,CB)*Dc(CB)+kG(CB,~CB)*Dx;
% R(CB)=Rx


