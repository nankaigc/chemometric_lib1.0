%% demo standardization by pca scores for corn protein
% written by Cheng Guo, College of Chemistry, Nankai University
clear;clc;close all;
%% load data
% parameter
data_name='corn';
idx_ana=3;
ratio=3;
method=2;
% load data
[ Data ]=prepare_data_latest(data_name,idx_ana,ratio,method );
X=[Data.Xm;Data.Xs;Data.Xs2];
[coeff,score,latent,tsquared,explained] = pca(X,'NumComponents',3);
%% score3
l=size(Data.Xm,1);
% 3D
figure;
scatter3(score(1:l,1),score(1:l,2),score(1:l,3),'b','filled');hold on;
scatter3(score(l+1:2*l,1),score(l+1:2*l,2),score(l+1:2*l,3),'r','filled');hold on;
scatter3(score(2*l+1:end,1),score(2*l+1:end,2),score(2*l+1:end,3),'g','filled');
ex1=eval(vpa(explained(1),4));ex2=eval(vpa(explained(2),3));ex3=eval(vpa(explained(3),2));
% xlabel({'1st Principal Component',ex1})
% ylabel({'2nd Principal Component',ex2})
% zlabel({'3rd Principal Component',ex3})
a=num2str(ex1);b=num2str(ex2);c=num2str(ex3);
xlabel({['PC1 ','(',a,'%',')']})
ylabel({['PC2 ','(',b,'%',')']})
zlabel({['PC3 ','(',c,'%',')']})
box on;
ConfidenceRegion(score(1:l,:));
ConfidenceRegion(score(l+1:2*l,:));
ConfidenceRegion(score(2*l+1:end,:));
legend1={'mp5';'mp6';'m5'};
legend(legend1,'Location','northwest');legend boxoff;
% 2D
figure;
scatter(score(1:l,1),score(1:l,2),'b','filled');hold on;
scatter(score(l+1:2*l,1),score(l+1:2*l,2),'r','filled');hold on;
scatter(score(2*l+1:end,1),score(2*l+1:end,2),'g','filled');
% xlabel({'1st Principal Component',ex1})
% ylabel({'2nd Principal Component',ex2})
xlabel({['PC1 ','(',a,'%',')']})
ylabel({['PC2 ','(',b,'%',')']})
box on;
ConfidenceRegion(score(1:l,1:2));
ConfidenceRegion(score(l+1:2*l,1:2));
ConfidenceRegion(score(2*l+1:end,1:2));
legend(legend1,'Location','northwest');legend boxoff;
clear a b c l ex1 ex2 ex3
%%
% subspace(score(1:54,:),score(109:end,:))
