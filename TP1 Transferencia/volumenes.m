%laploteadora

stepx=0.1;
stepy=0.2;

xvol=0:stepx:0.5;
yvol=0:stepy:1;

[Xv,Yv]=meshgrid(xvol,yvol);

hold on
plot(Xv,Yv,'k.')

line([0,0.5,0.5,0,0],[0,0,1,1,0])

xvol=xvol(1:end-1)+stepx/2;
yvol=yvol(1:end-1)+stepy/2;
for lx=xvol
    line([lx lx],[0 1])
end

for ly=yvol
    line([0 0.5],[ly ly])
end

xlim([-.05 0.55])
ylim([-.05 1.05])
daspect([1 1 1])
