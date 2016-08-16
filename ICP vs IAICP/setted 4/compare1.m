clear all
close all
clc
%%%%%%%%%%%%%%
% Loading files
load bdg1_tr_norm.txt
load bdg2_tr_norm.txt
% Zscore intensity normalization
Intensity_Data=bdg1_tr_norm(:,4);
Intensity_Model=bdg2_tr_norm(:,4);
ZnormInt_Data= ZscoreNormalizeIntensity(bdg1_tr_norm);
ZnormInt_Model= ZscoreNormalizeIntensity(bdg2_tr_norm);

% Run Intensity Augmented ICP
clear bdg1_tr_norm bdg2_tr_norm
D=ZnormInt_Data';
M=ZnormInt_Model';
m1=length(M);
d1=length(D);
% Run ICP (standard settings)
iter=50;
%tr1=[mode(D(1,:))-mode(M(1,:));mode(D(2,:))-mode(M(2,:));mode(D(3,:))-mode(M(3,:));0];
tr1=[2;10;2;0];
inT=repmat(tr1,1,d1);
D=D-inT;

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
dlmwrite('rgd_wt_IAICP.txt',registered,'newline','pc');
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
iter=50;
% tr1=[mode(D(1,:))-mode(M(1,:));mode(D(2,:))-mode(M(2,:));mode(D(3,:))-mode(M(3,:))];
tr1=[2;10;2];
inT=repmat(tr1,1,d1);
D=D-inT;

[Ricp Ticp ER2] = PointICP(M, D, iter,'Matching','kDtree');
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
dlmwrite('rgd_wt_ICP.txt',registered,'newline','pc');
% RMS error
fprintf('RMS error= %f',ER2);

ax=1:length(ER11);
pk=plot(ax,ER11,'--');
set(pk,'Color','red','LineWidth',1.25)
hold on
qk=plot(ax,ER21);
set(qk,'Color','blue','LineWidth',1.25)
xlabel('Iterations ');
ylabel('RMS Error');
legend('Intensity Augmented ICP','Geometric ICP');



