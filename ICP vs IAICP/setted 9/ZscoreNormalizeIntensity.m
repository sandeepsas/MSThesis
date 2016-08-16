function [DataOutput]= ZscoreNormalizeIntensity(DataInput)

A=DataInput(:,4);
mean_i=mean(A);
std_i=std(A);
N1=A-mean_i;
N2=N1/std_i;
DataOutput(:,4)=N2;
DataOutput(:,1:3)=DataInput(:,1:3);
