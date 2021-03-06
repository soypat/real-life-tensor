function [R] = qapply(q,x,nodes,nodestart,nodeend,R,ndof)
w=nodes(nodeend,1)-nodes(nodestart,1);
h=nodes(nodeend,2)-nodes(nodestart,2);
L=sqrt(w^2 +h^2);
theta=atan(h/w);
q1=q(x);
q2=q(x+L);
s=sin(theta);
c=cos(theta);
% q1x=s*q1;
% q2x=c*q2;

DeltaQ=q2-q1;
Atri=L/2*DeltaQ;
Acuadr=L*q1;
uposition1=nodestart*ndof-(ndof-1);
uposition2=nodeend*ndof-(ndof-1);
mposition1=uposition1+ndof-1;
mposition2=uposition2+ndof-1;
momsign=-sign(w*c*q1+s*q1*h);
if sign(DeltaQ)>=0
    R(uposition1)=R(uposition1)+s*(Acuadr/2+Atri/3);
    R(uposition2)=R(uposition2)+s*(Acuadr/2+Atri*2/3);
    if ndof>2
        R(uposition1+1)=R(uposition1+1)+c*(Acuadr/2+Atri/3);
        R(uposition2+1)=R(uposition2+1)+c*(Acuadr/2+Atri*2/3);
    end
    Mtri=L^2/12*abs(DeltaQ);
    Mcuadr=L^2*abs(q1)/24;
    R(mposition1)=R(mposition1)+momsign*(Mcuadr+.4*Mtri);
    R(mposition2)=R(mposition2)-momsign*(Mcuadr+.6*Mtri);
    return
else
    R(uposition1)=R(uposition1)+s*(Acuadr/2+Atri*2/3);
    R(uposition2)=R(uposition2)+s*(Acuadr/2+Atri/3);
    if ndof>2
        R(uposition1+1)=R(uposition1+1)+c*(Acuadr/2+Atri*2/3);
        R(uposition2+1)=R(uposition2+1)+c*(Acuadr/2+Atri/3);
    end
    Mtri=L^2/12*abs(DeltaQ);
    Mcuadr=L^2*abs(q2)/24;
    R(mposition1)=R(mposition1)+momsign*(Mcuadr+.6*Mtri);
    R(mposition2)=R(mposition2)-momsign*(Mcuadr+.4*Mtri);
        return
end 
end

