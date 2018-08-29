function [T] = Tvuw(u,v)
%TVUW This is easy
    L=zeros(3);
    w=cross(u,v);
    n=norm(u);
%     L([1 1 1],)
    L(:,1)=u/norm(u);
%     L(1,1)=u(1)/n;
%     L(2,1)=u(2)/n;
%     L(3,1)=u(3)/n;
%     n=norm(v);
    L(:,2)=v/norm(v);
    L(:,3)=w/norm(w);
%     L(1,3)=w(1)/n;
%     L(2,3)=w(2)/n;
%     L(3,3)=w(3)/n;
    T=zeros(12);
    T([1:3],[1:3])=L;
    T([4:6],[4:6])=L;
    T([7:9],[7:9])=L;
    T([10:12],[10:12])=L;
end

