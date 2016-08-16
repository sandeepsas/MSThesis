clear all
close all
clc
%%%%%%%%%%%%%%
% Loading files
load w_scan1.txt
load w_scan2.txt
% Zscore intensity normalization
Intensity_Data=w_scan1(:,4);
Intensity_Model=w_scan2(:,4);
ZnormInt_Data= ZscoreNormalizeIntensity(w_scan1);
ZnormInt_Model= ZscoreNormalizeIntensity(w_scan2);

% Run Intensity Augmented ICP
clear w_scan1 w_scan2
D=ZnormInt_Data';
M=ZnormInt_Model';
m1=length(M);
d1=length(D);
% Run ICP (standard settings)
iter=50;
% tr1=[mean(D(1,:))-mean(M(1,:));mean(D(2,:))-mean(M(2,:));mean(D(3,:))-mean(M(3,:));0];
% inT=repmat(tr1,1,d1);
% D=D-inT;

[Ricp Ticp ER1] = IAPICP(M, D, iter,0);
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
% tr1=[mean(D(1,:))-mean(M(1,:));mean(D(2,:))-mean(M(2,:));mean(D(3,:))-mean(M(3,:));0];
% tr1=[-0.1;10;0.1;0];
% inT=repmat(tr1,1,d1);
% D=D-inT;

[Ricp Ticp ER2] = PointICP(M, D, iter,'Matching','kDtree','Minimize','plane');
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


ER11=ER1(1,2:end);
ER21=ER2(2:end,:)';

ax=1:length(ER11);
pk=plot(ax,ER11,'--');
set(pk,'Color','red','LineWidth',1.25)
hold on
qk=plot(ax,ER21);
set(qk,'Color','blue','LineWidth',1.25)
xlabel('Iterations ');
ylabel('RMS Error');
legend('Intensity Augmented ICP','Geometric ICP');



