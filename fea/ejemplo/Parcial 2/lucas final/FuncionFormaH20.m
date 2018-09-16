function [N D]=FuncionFormaH20(ksi,eta, zeta)
syms r s t ksi eta zeta
X=[ 1, r, s, t, r^2, s^2, t^2, r*s, s*t, r*t, r^2*s, r^2*t,s^2*r,s^2*t,t^2*r,t^2*s,r*s*t,r^2*s*t,s^2*r*t,t^2*r*s];
Nodos=[-1    -1   -1
        1    -1   -1
        1     1   -1
       -1     1   -1
       -1    -1    1
        1    -1    1
        1     1    1
       -1     1    1
        0    -1   -1
        1     0   -1
        0     1   -1
       -1     0   -1
       -1    -1    0
        1    -1    0
        1     1    0
       -1     1    0
        0    -1    1
        1     0    1
        0     1    1
       -1     0    1];
   A=sym(zeros(length(X)));
   for i=1:length(X)
       A(i,:)=subs(X,[r s t],[Nodos(i,1) Nodos(i,2) Nodos(i,3)]);
   end
A=double(A);
N=X/(A);
% Check=zeros(length(X));
%    for i=1:length(X)
%        Check(i,:)=subs(N,[r s t],[Nodos(i,1) Nodos(i,2) Nodos(i,3)]);
%    end
%    Check;
   D=sym(zeros(2,length(X)));
   D(1,:)=diff(N,r);
   D(2,:)=diff(N,s);
   D(3,:)=diff(N,t);
   N=double(subs(N,[r s t],[ksi eta zeta]));
   D=double(subs(D,[r s t],[ksi eta zeta]));
end