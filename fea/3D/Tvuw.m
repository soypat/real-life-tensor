function [T] = Tvuw(u,v)
%TVUW This is easy
    L=zeros(3);
    w=cross(u,v);
    n=norm(u);
    L(:,1)=u/norm(u);
    L(:,2)=v/norm(v);
    L(:,3)=w/norm(w);
    T=zeros(12);
    T=blkdiag(L,L,L,L);
end

