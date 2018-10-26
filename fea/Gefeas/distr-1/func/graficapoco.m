function [] = graficapoco(nod,elenod,eletype,Ie)
clf
Imax=max(Ie);
Imin=min(Ie);
if Imax==Imin
    Imax=Imax+100;
    Imin=Imin-20;
end
thickness=@(i) (4/(Imax-Imin))*i+(1-4*Imin/(Imax-Imin));
[Ne,~]=size(elenod);
hinges=zeros(sum(ismember(eletype,3))+sum(ismember(eletype,4)),2);
hingecount=1;
for i = 1:Ne
    hold on
    xv=[nod(elenod(i,1),1) nod(elenod(i,2),1)];
    yv=[nod(elenod(i,1),2) nod(elenod(i,2),2)];
    myline=line(xv,yv);
    switch eletype(i)
        case 1
            color='k';
        case 2
            color=[0.5 0.5 0.5];
        case 3
            color=[0.5 0.5 0.5];
            hinges(hingecount,[1 2])=[nod(elenod(i,1),1),nod(elenod(i,1),2)];
            hingecount=hingecount+1;
        case 4
            color=[0.5 0.5 0.5];
            hinges(hingecount,[1 2])=[nod(elenod(i,2),1),nod(elenod(i,2),2)];
            hingecount=hingecount+1;
        case 11
            color=[0 45 119]/256;
        case 22
            color=[0 119 45]/256;
        case 33
            color=[0 119 45]/256;
            hinges(hingecount,[1 2])=[nod(elenod(i,1),1),nod(elenod(i,1),2)];
            hingecount=hingecount+1;
        case 44
            color=[0 119 45]/256;
            hinges(hingecount,[1 2])=[nod(elenod(i,2),1),nod(elenod(i,2),2)];
            hingecount=hingecount+1;
    end
    set(myline,'LineWidth',thickness(Ie(i)),'Color',color)
    
end
scatter(nod(:,1),nod(:,2),'.k')
scatter(hinges(:,1),hinges(:,2),'ok')
end