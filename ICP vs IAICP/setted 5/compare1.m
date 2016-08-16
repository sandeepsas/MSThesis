clear all
close all
clc
%%%%%%%%%%%%%%
% Loading files
load scan3_p2_tr.txt
load scan3_p2.txt
% Zscore intensity normalization
Intensity_Data=scan3_p2_tr(:,4);
Intensity_Model=scan3_p2(:,4);
 ZnormInt_Data= ZscoreNormalizeIntensity(scan3_p2_tr);
 ZnormInt_Model= ZscoreNormalizeIntensity(scan3_p2);
% Scaled_Intensity_Data=scan3_p2_tr(:,4)/20;
% Scaled_Intensity_Model=scan3_p2(:,4)/20;
% % 
% ZnormInt_Data= [scan3_p2(:,1:3),Scaled_Intensity_Data];
% ZnormInt_Model=[scan4_p2(:,1:3),Scaled_Intensity_Model];

% Run Intensity Augmented ICP
clear scan3_p2 scan4_p2
D=ZnormInt_Data';
M=ZnormInt_Model';
m1=length(M);
d1=length(D);
% Run ICP (standard settings)
iter=30;
% tr1=[mean(D(1,:))-mean(M(1,:));mean(D(2,:))-mean(M(2,:));mean(D(3,:))-mean(M(3,:));0];
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


D=ZnormInt_Data(:,1:3)';
M=ZnormInt_Model(:,1:3)';
m1=length(M);
d1=length(D);
% Run ICP (standard settings)
iter=30;
% tr1=[mean(D(1,:))-mean(M(1,:));mean(D(2,:))-mean(M(2,:));mean(D(3,:))-mean(M(3,:));0];
% tr1=[-0.1;10;0.1;0];
% inT=repmat(tr1,1,d1);
% D=D-inT;

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



