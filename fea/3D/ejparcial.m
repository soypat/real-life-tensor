%parcial
nodaux=linspace(0,1,10);
nod=[0 0 0;.1 0 0;.2 0 0;.3 0 0;.4 0 0;.5 0 0;.6 0 0;.7 0 0;.8 0 0;.9 0 0;1 0 0;1 10 0];
nod=nod*1000;
enod=[1 2;2 3;3 4;4 5;5 6;6 7;7 8;8 9;9 10;10 11;11 12];

fnod=@(n) [n*6-5 n*6-4 n*6-3 n*6-2 n*6-1 n*6];
[Ne,~]=size(enod);
[N,~]=size(nod);

Le=zeros(Ne,1);

Ts={};
kG=zeros(6*N);
A=7.64e2;
L=1000;
E=210e6;
K=0.721e4;
Iz=80.1e4;
Iy=8.49e4;


for i=1:Ne
    ns=enod(i,1);
    ne=enod(i,2);
    lx=nod(ne,1)-nod(ns,1);
    ly=nod(ne,2)-nod(ns,2);
    lz=nod(ne,3)-nod(ns,3);
    Le(i)=sqrt(lx^2+ly^2+lz^2);%Long.
    px=[lx ly lz]/Le(i);
    vy=[0 1 pi];
    pz=cross(px,vy)/norm(cross(px,vy));
    py=cross(pz,px);
    C=[px;py;pz;];%Para barras no importa orientacion azimuth (acimutal)
    T=blkdiag(C,C,C,C);
    Ts=[Ts T];
%     klocal=Kvuw(Ee(i),Ae(i),0,0,0,0,Le(i));
    klocal=Kvuw(E,A,Iz,Iy,K,.3,Le(i));
    krotada=T'*klocal*T;
    kG([fnod(enod(i,1)) fnod(enod(i,2))],[fnod(enod(i,1)) fnod(enod(i,2))])=...
      kG([fnod(enod(i,1)) fnod(enod(i,2))],[fnod(enod(i,1)) fnod(enod(i,2))])+krotada;
end

CB=false(N*6,1);
CB(fnod(1))=true;
CB(fnod(12))=true;
Kr=kG(~CB,~CB);
R=zeros(6*N,1);
q=-2000;

for i=1:10
    R(i*6-2)=q*Le(2)/2+R(i*6-2);
    R((i+1)*6-2)=q*Le(i)/2+R((i+1)*6-2);
end
% R((i+1)*6-2)=q*Le(2)/2;
% R(6*10-2)=-200000;
% 
R=zeros(6*N,1);

R(11*6-2)=-2000*Le(2)*(Ne-1)/2;
R(1*6-2)=-2000*Le(2)*(Ne-1)/2;
F=R(~CB);

U=Kr\F;
D=zeros(6*N,1);
D(~CB)=U;
rads=D(fnod(11));
rads=rads(4);
disp(rads)
disp(D(fnod(12)))