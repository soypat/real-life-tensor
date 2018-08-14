function [R] = qapplyx(q,x,L,nodestart,nodeend,R,ndof)
q1=q(x);
q2=q(x+L);
DeltaQ=q2-q1;
Atri=L/2*DeltaQ;
Acuadr=L*q1;
uposition1=nodestart*ndof-(ndof-1);
uposition2=nodeend*ndof-(ndof-1);
mposition1=uposition1+ndof-1;
mposition2=uposition2+ndof-1;
if sign(DeltaQ)>=0
        R(uposition1)=R(uposition1)+Acuadr/2+Atri/3;
        R(uposition2)=R(uposition2)+Acuadr/2+Atri*2/3;
        Mtri=L^2/12*DeltaQ;
        Mcuadr=L^2*q1/24;
        R(mposition1)=R(mposition1)+Mcuadr+.4*Mtri;
        R(mposition2)=R(mposition2)-Mcuadr-.6*Mtri;
        return
else
        R(uposition1)=R(uposition1)+Acuadr/2+Atri*2/3;
        R(uposition2)=R(uposition2)+Acuadr/2+Atri/3;
        Mtri=L^2/12*abs(DeltaQ);
        Mcuadr=L^2*q2/24;
        R(mposition1)=R(mposition1)+Mcuadr+.6*Mtri;
        R(mposition2)=R(mposition2)-Mcuadr-.4*Mtri;
        return
end 
end

