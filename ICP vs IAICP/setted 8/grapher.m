clc
clear all
close all

ER11=[0.3563,0.8976,0.0326,0.0716];
ER21=[0.0577,0.8256,0.0261,0.0318];

ax=1:length(ER11);
pk=plot(ax,ER11,'--x');
set(pk,'Color','blue','LineWidth',1.25)
hold on
qk=plot(ax,ER21,'--x');
set(qk,'Color','red','LineWidth',1.25)
xlabel('Data Set ');
ylabel('RMS Error');
legend('Geometric ICP','Intensity Augmented ICP');