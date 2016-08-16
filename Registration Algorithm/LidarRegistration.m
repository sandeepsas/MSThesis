% Clearing workspace
clc
clear all
close all
% Loading files
load Scan3_slave.txt
load Scan3_master.txt
% Zscore intensity normalization
Scan3_slave(:,4)=Scan3_slave(:,4)+randi(10,size(Scan3_slave,1),1);
Scan3_master(:,4)=Scan3_master(:,4)+randi(10,size(Scan3_master,1),1);
Intensity_Data=Scan3_slave(:,4);
Intensity_Model=Scan3_master(:,4);
ZnormInt_Data= ZscoreNormalizeIntensity(Scan3_slave);
ZnormInt_Model= ZscoreNormalizeIntensity(Scan3_master);

% Run Intensity Augmented ICP
clear Scan3_slave Scan3_master
D=ZnormInt_Data';
M=ZnormInt_Model';
m1=length(M);
d1=length(D);
% Run ICP (standard settings)
iter=30;
tr1=[mean(D(1,:))-mean(M(1,:));mean(D(2,:))-mean(M(2,:));mean(D(3,:))-mean(M(3,:));0];
tr1=[-0.1;10;0.1;0];
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
dlmwrite('coarsereg.txt',registered,'newline','pc');
% RMS error
fprintf('RMS error= %f',ER1);

% Geometric ICP

data=Dicp;
model=M;
Intensity_Data=data(:,4);
Intensity_Model=model(:,4);
Intensity=[Intensity_Model;Intensity_Data];
M=model(:,1:3);
D=data(:,1:3);
M=M';D=D';
clear model data
[T R Ricp1 Ticp1 ER2 t ] = GeometricICP(M, D, 20,0);
n=length(D);
Dicp = Ricp1 * D + repmat(Ticp1, 1, n);

rg=horzcat(M,Dicp);
rg=rg';
rg(:,4)=Intensity;
dlmwrite('finereg.txt',rg,'newline','pc');
% dlmwrite('model_bdg_plane.txt',M,'newline','pc');
% dlmwrite('data_bdg_plane.txt',Dicp,'newline','pc');
Rtot=(Ricp*Ricp1);
Ttot=-1*(Ricp1*Ticp+Ticp1);
Rtot=inv(Rtot);


ER=[ER1,ER2];
ax=1:length(ER);
plot(ax,ER,'--x');
xlabel('Iterations ');
ylabel('RMS Error');
legend('Intensity Augmented Registration');

