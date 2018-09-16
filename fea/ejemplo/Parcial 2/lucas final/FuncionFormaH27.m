function [N D]=FuncionFormaH20(ksi,eta, zeta)
syms x y z ksi eta zeta
X=[ 1, x, y, z, x^2, y^2, z^2, x*y, y*z, x*z, x^2*y, x^2*z,y^2*x,y^2*z,z^2*x,z^2*y,x*y*z,x^2*y*z,y^2*x*z,z^2*x*y,x^2*y^2,y^2*z^2,z^2*x^2,x^2*y^2*z,y^2*z^2*x,z^2*x^2*y,x^2*y^2*z^2];
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
       -1     0    1
        0     0   -1
       -1     0    0
        0    -1    0
        1     0    0
        0     1    0
        0     0    1
        0     0    0];
   A=sym(zeros(length(X)));
   for i=1:length(X)
       A(i,:)=subs(X,[x y z],[Nodos(i,1) Nodos(i,2) Nodos(i,3)]);
   end
A=double(A);
N=X/(A);
% Check=zeros(length(X));
%    for i=1:length(X)
%        Check(i,:)=subs(N,[x y z],[Nodos(i,1) Nodos(i,2) Nodos(i,3)]);
%    end
%    Check;
   D=sym(zeros(2,length(X)));
   D(1,:)=diff(N,x);
   D(2,:)=diff(N,y);
   D(3,:)=diff(N,z);
   N=double(subs(N,[x y z],[ksi eta zeta]));
   D=double(subs(D,[x y z],[ksi eta zeta]));
end