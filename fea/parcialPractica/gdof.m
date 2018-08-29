function [outputArg1,outputArg2] = gdof(nod,elenod, eletype)
[Ne,~]=size(elenod);
[N,~]=size(nod);
hinges=length(eletype(eletype>2));
Le=zeros(Ne,1);
phide=zeros(Ne,1);
Ndof=N*3+hinges;
hingedlist=zeros(1,hinges);
j=0;
Nt=N+hinges;
nodeDofs=zeros(Nt,3);
nudos=zeros(Nt,1);
for i=1:N
    nodeDofs(i,:)=[i*3-2 i*3-1 i*3];
end
for i = 1:Ne
    nodestart=elenod(i,1);
    nodeend=elenod(i,2);
    lx=nod(nodeend,1)-nod(nodestart,1);
    ly=nod(nodeend,2)-nod(nodestart,2);
    Le(i)=sqrt(lx^2+ly^2);
    phide(i)=atan2d(ly,lx);
    
end

