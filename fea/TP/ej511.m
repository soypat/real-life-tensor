nodes=[0 0;
    40 0;
    70 40;
    20 40];
elenod=[1 4;
        2 4;
        3 4];
[Ne,~]=size(elenod);
elementos=genelementos(elenod);
Le=zeros(Ne,1);
ndof=3;

for i = 1:Ne
    nodestart=elenod(i,1);
    nodeend=elenod(i,2);
    lx=nodes(nodeend,1)-nodes(nodestart,1);
    ly=nodes(nodeend,2)-nodes(nodestart,2);
    Le(i)=sqrt(lx^2+ly^2);
end