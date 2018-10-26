%empiezo 10:47

nod=[0 0 0;-500  0 -866.03;500 0 -866.03;-500 1500 -866.03;%nod4
    0 1500 0;500 1500 -866.3;0 1500+816.5 -577.35];
enod=[1 5;2 4;3 6;4 6;4 5;6 5;4 7;6 7;5 7];
Ne=size(enod,1);
N=size(nod,1);
scatter3(nod(:,1),nod(:,2),nod(:,3))
% for i=1:Ne
%     hold on
%     plot3()
% end
Ts={};
kcell={};
Le=zeros(Ne,1);
b=20;
h=150;
d=50;
fnod=@(n) [n*6-5 n*6-4 n*6-3 n*6-2 n*6-1 n*6];
fdof=@(n,i) n*6-5+i;
kG=zeros(N*6);
for i=1:Ne
    ns=enod(i,1);
    ne=enod(i,2);
    lx=nod(ne,1)-nod(ns,1);
    ly=nod(ne,2)-nod(ns,2);
    lz=nod(ne,3)-nod(ns,3);
    Le(i)=sqrt(lx^2+ly^2+lz^2);
    px=[lx ly lz]/Le(i);
    switch i
        case 1
            vy=[0 0 1];
        case {2,3}
            vy=[1 0 0];
        otherwise
            vy=[0 pi 0];
    end
    nu=.3;
    switch i
        case {1,2,3}
            Iz=b*h^3/12;
            Iy=h*b^3/12;
            A=b*h;
            K=(Iz+Iy)/2;
        otherwise
            Iz=pi*d^4/64;
            Iy=Iz;
            A=pi*d^2/4;
            K=Iz*2;
    end
    pz=cross(px,vy)/norm(cross(px,vy));
    py=cross(pz,px);
    C=[px' py' pz'];
    T=blkdiag(C,C,C,C);
    Ts=[Ts T];
    klocal=Kvuw(210e3,A,Iz,Iy,K,nu,Le(i));%Mpa,mm,mm^2,mm^4
    kcell=[kcell klocal];
    krot=T'*klocal*T;
    index=[fnod(ns) fnod(ne)];
    kG(index,index)=kG(index,index)+krot;
end

CB=false(6*N,1);
CB([fnod(1) fnod(2) fnod(3)])=true;
A=b*h;
Dt=110;
E=210e3;
alphas=1.2e-5;
Fterm=A*E*alphas*Dt;
 
R=zeros(N*6,1);
R(fnod(7))=[1.5e3 0 0 0 0 0];
R(fnod(6))=[0 Fterm 0 0 0 0];
sigt=Fterm/A;
Kr=kG(~CB,~CB);
F=R(~CB);

U=Kr\F;
D=zeros(N*6,1);
D(~CB)=U;
U6=D(fnod(6));
U3=D(fnod(3));
ulocal3=[U3; U6];
Udt3=U3(2);
Udt6=U6(2);
d=[fdof(3,2) fdof(6,2)]
sig3=E*(Udt6-Udt3)/Le(3)+sigt;
F3=sig3*A;


%          u         v                 w       t1        t2        t3     
Bz=@(x,L) [-1/L -6/L^2+12*x/L^3 0 0   0 -4/L+6*x/L^2 ...
           1/L  6/L^2-12*x/L^3  0 0   0 -2/L+6*x/L^2];
By=@(x,L) [-1/L 0 -6/L^2+12*x/L^3 0   -4/L+6*x/L^2 0 ...
           1/L 0 6/L^2-12*x/L^3   0   -2/L+6*x/L^2 0 ];
kurv=Bz(1,1.5)*ulocal3;
kurv*E*Iz