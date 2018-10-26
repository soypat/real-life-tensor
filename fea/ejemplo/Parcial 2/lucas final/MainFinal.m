
for i=1:11
    tic
    i
    Espesor1=0.2*i-0.2;
   [avgStress8,AvgVonMisses,Dnodes,nodes,elements]=FuncFinal2(Espesor1);
    Avg(i).Stress8=avgStress8;
    Avg(i).VonMisses=AvgVonMisses;
    D(i).Dnodes=Dnodes;
   
    toc
end
for j=1:11
AvgStress8=zeros(38509,6);
for i=1:8400
    AvgStress8(elements(i,:),1)=Avg(j).Stress8(i,:,1);
    AvgStress8(elements(i,:),2)=Avg(j).Stress8(i,:,2);
    AvgStress8(elements(i,:),3)=Avg(j).Stress8(i,:,3);
    AvgStress8(elements(i,:),4)=Avg(j).Stress8(i,:,4);
    AvgStress8(elements(i,:),5)=Avg(j).Stress8(i,:,5);
    AvgStress8(elements(i,:),6)=Avg(j).Stress8(i,:,6);
end
Avg(j).AvgStress8=AvgStress8;
end
for j=1:11
AvgStress8=zeros(38509,6);
for i=1:8400
    AvgStress8(elements(i,:),1)=Avg(j).Stress8(i,:,1);
    AvgStress8(elements(i,:),2)=Avg(j).Stress8(i,:,2);
    AvgStress8(elements(i,:),3)=Avg(j).Stress8(i,:,3);
    AvgStress8(elements(i,:),4)=Avg(j).Stress8(i,:,4);
    AvgStress8(elements(i,:),5)=Avg(j).Stress8(i,:,5);
    AvgStress8(elements(i,:),6)=Avg(j).Stress8(i,:,6);
end
Avg(j).AvgStress8=AvgStress8;
end
Name='AceroCobreTermico';
    mkdir(Name);
   filename=sprintf('%s',strcat('C:\Users\Lucas\Desktop\Elementos Finitos\Final\',Name,'\AceroCobre2.mat'));
    save(filename,'Avg','D','nodes','elements');
   