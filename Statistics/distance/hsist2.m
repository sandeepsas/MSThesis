clc

close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mean1=44;
mean2=53;
sig1=8;
sig2=9;


num1=2*sig1*sig2;
den1=sig1^2+sig2^2;
unrt=sqrt(num1/den1);

num2=(mean1-mean2)^2;
e_up=-0.25*(num2/den1);
sub1=e_up*unrt;
dist2=sqrt(1-sub1)
