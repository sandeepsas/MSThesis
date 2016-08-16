% Clearing workspace
clc
clear all
close all
% Loading files
load bdg2_slave.txt
load bdg2_master.txt
% Zscore intensity normalization
Intensity_Data=bdg2_slave(:,4);
Intensity_Model=bdg2_master(:,4);
ZnormInt_Data= ZscoreNormalizeIntensity(bdg2_slave);
ZnormInt_Model= ZscoreNormalizeIntensity(bdg2_master);

% Run Intensity Augmented ICP
clear bdg2_slave bdg2_master
D=ZnormInt_Data';
M=ZnormInt_Model';
m1=length(M);
d1=length(D);
% Run ICP (standard settings)
iter=30;
% tr1=[mean(D(1,:))-mean(M(1,:));mean(D(2,:))-mean(M(2,:));mean(D(3,:))-mean(M(3,:));0];
% tr1=[-0.1;10;0.1;0];
% inT=repmat(tr1,1,d1);
% D=D-inT;

[Ricp Ticp ER1] = IntensityAugmentedICP(M, D, iter,0);
%
% Transform data-matrix using ICP result
D=D(1:3,:);
% M=M(1:3,:);
Dicp = Ricp * D + repmat(Ticp, 1, d1);
M=M';Dicp=Dicp';
M(:,4)=Intensity_Model;
Dicp(:,4)=Intensity_Data;
% dlmwrite('model_bdg.txt',M,'newline','pc');
% dlmwrite('data_bdg.txt',Dicp,'newline','pc');
registered=vertcat(M,Dicp);
%registered=registered';
dlmwrite('rgd_bdg_IAICP.txt',registered,'newline','pc');
% RMS error
fprintf('RMS error= %f',ER1);

ER=ER1;
ax=1:length(ER);
plot(ax,ER,'--x');
xlabel('Iterations ');
ylabel('RMS Error');
legend('Intensity Augmented Registration');

