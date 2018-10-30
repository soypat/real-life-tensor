%% 
T=5000;
c1=1;
c2=1;
lam1=3e-5;
lam2=3e-4;
 planck=@(L) c1./(L.^5.*(exp(c2./L./T)-1));
 n=16;
 lam=linspace(lam1,lam2,n);
 
 vo=integral(planck,lam1,lam2);
 done=0;
 while done==0
     n=2*n;
     lam=linspace(lam1,lam2,n);
     vn=integral(planck,lam1,lam2);
     if abs((vn-vo)/vo)<1e-3
         done=1;
         vf=vn;
     end
 end