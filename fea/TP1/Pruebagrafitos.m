clear
clc
% clf
close all

fl1=[10000,-1200,750,-10000,1200,800];
 %fl1=[1200,750,-1200,800];
L=3;
y=0.1;
b=0.2;
A=b*2*y;
I=(b*(2*y)^3)/12;
l=length(0:L);
E=20e9;
% axis tight

%GrafitodoF(L,L/3,2*L,fl1,-fl1,2*fl1)
St=Calcutodo(A,I,y,b,E,fl1);