P=3200;
E=30e6;
h_nom=10;%in
b_nom=10;
Fa=-P*15/16;
Fb=-Fa;
nodes=[0 0;
    3 0;
    3 3];

elementos=[1 2 3 4 5 6;
    1 2 3 7 8 9];

[N,placeholder]=size(nodes);
Ne=5;
ndof=3;
Ndof=N*ndof;
barras=[2];
kG=zeros(Ndof);
I1=.0004;
A1=0.02;
A3=0.001;
A4=0.003;
A5=0.001;



for i  = 1:Ne
    nodestart=elementos(i,ndof)/ndof;
    nodeend=elementos(i,ndof*2)/ndof;
    lx=nodes(nodeend,1)-nodes(nodestart,1);
    ly=nodes(nodeend,2)-nodes(nodestart,2);
    L=sqrt(lx^2+ly^2);
    phi=atan(ly/lx);
    phid=atand(ly/lx);
    if sum(ismember(barras,i))>0
         A=1e-3;
         E=210e9;
         klocal=Kb(E,A,L,phid);
         kG(elementos(i,:),elementos(i,:))=kG(elementos(i,:),elementos(i,:))+klocal;
    else
        A=2e-3;
        I=5e-5;
        klocal=Kv(30e9,A1,I1,L,phid);
        kG(elementos(i,:),elementos(i,:))=kG(elementos(i,:),elementos(i,:))+klocal;
    end
end


