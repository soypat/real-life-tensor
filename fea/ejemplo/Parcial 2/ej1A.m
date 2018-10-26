clc
clear

syms x y real

pol = [ 1 x y x^2 x*y];
A = [ 1 -1 -1 1 1 
      1 1 -1 1 -1
      1 1 1 1 1
      1 -1 1 1 -1
      1 0 0 0 0];
 N = pol / (A);
 dNx = diff(N,x);

 dNy = diff(N,y);
 display(N)
 display(dNx)
 display(dNy)
 sum(N)
 sum(dNx)
 sum(dNy)