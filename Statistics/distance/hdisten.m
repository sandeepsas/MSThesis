clc

close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mean1=input('mean1','s');
mean2=input('mean2','s');
sig1=input('sigma1','s');
sig2=input('sigma2','s');

 mean1=str2num(mean1);
 mean2=str2num(mean2);
 sig1=str2num(sig1);
 sig2=str2num(sig2);

num1=2*sig1*sig2;
den1=sig1^2+sig2^2;
unrt=sqrt(num1/den1);

num2=(mean1-mean2)^2;
e_up=-0.25*(num2/den1);
sub1=e_up*unrt;
dist1=sqrt(1-sub1)