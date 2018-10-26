syms x y b h real
Nvec=[1-x/b-y/h,0, y/h,0, x/b];

Ncst=[Nvec,0;0,Nvec];
% Bcst=doper(N);

%% Nlst?
nod=[0 0;b 0;0 h;b/2 h/2;0 h/2;b/2 0];
syms a1 a2 a3 a4 a5 a6
u=@(x,y) a1+a2*x+a3*y;
v=@(x,y) a4+a5*x+a6*y;
amat=[a1;a2;a3];
A=[1 0 0;1 b 0;1 0 h];
U=[1;0;0];
au1=A\U
A*au1
%% constitutiva
nu=.2;
E=1;
C = E/(1 - nu^2)*[ 1.0      nu        0.0
                   nu       1.0       0.0
                   0.0      0.0    (1 - nu)/2 ];%
%

