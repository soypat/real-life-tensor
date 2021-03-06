%% Rotador de vigas 3D
function Ke=rotar(Kl,p1,p2,p3)
% p1=[0 0 0];
% p2=[1 1 0];
% p3=[0 0 1];
v1 = p2-p1;
vd1 = v1/norm(v1);
vp = p3-p1;
vd3 = cross(vd1,vp)/norm(cross(vd1,vp));
vd2 = cross(vd3,vd1);
lambda = [ vd1 ; vd2 ; vd3 ];
T=zeros(12);
for i = 1:4
    loc = 3*(i-1)+(1:1:3);
    T(loc,loc) = lambda;
end
Ke = T.'*Kl*T;
end
% a=lambda*v1';