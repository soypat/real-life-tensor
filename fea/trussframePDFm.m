P=3200;
E=30e6;
h_nom=10;%in
b_nom=10;
Fa=-P*15/16;
Fb=-Fa;
nodes=[-4 3;
    0 0;
    0 3;
    4 3];

elementos=[1 2 3 7 8 9;%viga1
    1 2 3 4 5 6;%barra 3
    4 5 6 7 8 9;%barra 4 UP
    7 8 9 10 11 12;% viga2
    4 5 6 10 11 12];%barra5
[N,placeholder]=size(nodes);
Ne=5;
ndof=3;
Ndof=N*ndof;
barras=[3 4 5];
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
         A=eval(sprintf('A%0.0f',i));
         E=200e9;
         klocal=Kb(E,A,L,phid);
         kG(elementos(i,:),elementos(i,:))=kG(elementos(i,:),elementos(i,:))+klocal;
    else
        klocal=Kv(30e9,A1,I1,L,phid);
        kG(elementos(i,:),elementos(i,:))=kG(elementos(i,:),elementos(i,:))+klocal;
    end
end


