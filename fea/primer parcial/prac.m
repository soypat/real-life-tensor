% practica
nod=[0 0 0;1 0 1;2 0 0;0 1 0;1 1 1;2 1 0];
elenod=[1 4;2 5;3 6;4 2;5 3;6 1;4 5;5 6;6 4];
Ne=size(elenod,1);
N=size(nod,1);
eletype=ones(1,Ne);
eletype([1 2 3])=[2 2 2];
px=cell(Ne,1);
py=cell(Ne,1);
pz=cell(Ne,1);
T=cell(Ne,1);

py{1}=[-1 0 0];
py{2}=[0 0 1];
py{3}=[-1 0 0];
Le=zeros(Ne,1);
kloc=cell(Ne,1);
E=200e9;
b=.01;
h=.08;
A=h*b;
Iz=b*h^3/12;
Iy=h*b^3/12;
K=A*(Iz+Iy)/((h+b)^2);
nu=.3;
fnod=@(n) [6*n-5 6*n-4 6*n-3 6*n-2 6*n-1 6*n];
fnodb=@(n) [6*n-5 6*n-4 6*n-3];
kG=zeros(N*6);
for i=1:Ne
    ns=elenod(i,1);
    ne=elenod(i,2);
    px{i}=(nod(ne,:)-nod(ns,:))/norm(nod(ne,:)-nod(ns,:));
    Le(i)=norm(nod(ne,:)-nod(ns,:));
    switch i
        case {1 2 3}
            pz{i}=cross(px{i},py{i})/norm(cross(px{i},py{i}));
            kloc{i}=Kvuw(E,A,Iz,Iy,K,nu);
            index=[fnod(ns) fnod(ne)];
            kG(index,index)=kG(index,index)+kloc{i};
        otherwise
            py{i}=[pi exp(1) 2]/norm([pi exp(1) 2]);
            pz{i}=cross(px{i},py{i})/norm(cross(px{i},py{i}));
            fullk=Kvuw(E,A,0,0,0,.3);
            kloc{i}=fullk([])
            index=[fnodb(ns) fnodb(ne)];
            kG(index,index)=kG(index,index)+kloc{i};
    end
    
    C=[px{i};py{i};pz{i}];
    T{i}=blkdiag(C,C,C,C);
    
end