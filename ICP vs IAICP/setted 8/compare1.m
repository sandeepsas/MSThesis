clear all
close all
clc
%%%%%%%%%%%%%%
% Loading files
load bdg2_slave.xyz
load bdg2_master.xyz
% Zscore intensity normalization
Intensity_Data=bdg2_slave(:,4);
Intensity_Model=bdg2_master(:,4);
ZnormInt_Data= ZscoreNormalizeIntensity(bdg2_slave(:,1:4));
ZnormInt_Model= ZscoreNormalizeIntensity(bdg2_master(:,1:4));

% Run Intensity Augmented ICP
clear bdg2_slave bdg2_master
D=ZnormInt_Data';
M=ZnormInt_Model';
m1=length(M);
d1=length(D);
% Run ICP (standard settings)
iter=30;
tr1=[median(D(1,:))-median(M(1,:));median(D(2,:))-median(M(2,:));median(D(3,:))-median(M(3,:));0];
inT=repmat(tr1,1,d1);
D=D-inT;

[Ricp Ticp ER1] = IntensityAugmentedICP(M, D, iter,0.4);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


D=ZnormInt_Data(:,1:3)';
M=ZnormInt_Model(:,1:3)';
m1=length(M);
d1=length(D);
% Run ICP (standard settings)
iter=30;
 tr1=[median(D(1,:))-median(M(1,:));median(D(2,:))-median(M(2,:));median(D(3,:))-median(M(3,:))];
% tr1=[-0.1;10;0.1;0];
 inT=repmat(tr1,1,d1);
 D=D-inT;

[Ricp Ticp ER2] = PointICP(M, D, iter,'Matching','kDtree','WorstRejection',0.4);
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
dlmwrite('rgd_bdg_ICP.txt',registered,'newline','pc');
% RMS error
fprintf('RMS error= %f',ER2);


ER11=ER1(1,2:end);
ER21=ER2(2:end,:)';

ax=1:length(ER11);
pk=plot(ax,ER11,'--x');
set(pk,'Color','red','LineWidth',1.25)
hold on
qk=plot(ax,ER21,'--x');
set(qk,'Color','blue','LineWidth',1.25)
xlabel('Iterations ');
ylabel('RMS Error');
legend('Intensity Augmented ICP','Geometric ICP');



