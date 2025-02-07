load('Revise_Exp3_ClusterEnhancement_kmodes.mat')
load('Revise_Exp3_ClusterEnhancement_CDCDR.mat')
load('Revise_Exp3_ClusterEnhancement_Het2Hom.mat')
load('Revise_Exp3_ClusterEnhancement_CMS.mat')
load('Exp1_CMIvsACC_NMI_ARI_Fscore.mat')

%% Enhancement_k-modes
Table_Enhancement_kmodes = [Enhancement_ACC_kmodes_mean Enhancement_NMI_kmodes_mean Enhancement_ARI_kmodes_mean Enhancement_FSC_kmodes_mean];
%% k-modes
Table_kmodes = [com_ACC_kmodes_mean com_NMI_kmodes_mean com_ARI_kmodes_mean com_FSC_kmodes_mean];
%% Enhancement_CDCDR
Table_Enhancement_CDCDR = [Enhancement_ACC_CDCDR_mean Enhancement_NMI_CDCDR_mean Enhancement_ARI_CDCDR_mean Enhancement_FSC_CDCDR_mean];
%% CDCDR
Table_CDCDR = [com_ACC_CDCDR_mean com_NMI_CDCDR_mean com_ARI_CDCDR_mean com_FSC_CDCDR_mean];
%% Enhancement_Het2Hom
Table_Enhancement_Het2Hom = [Enhancement_ACC_Het2Hom_mean Enhancement_NMI_Het2Hom_mean Enhancement_ARI_Het2Hom_mean Enhancement_FSC_Het2Hom_mean];
%% Het2Hom
Table_Het2Hom = [com_ACC_Het2Hom_mean com_NMI_Het2Hom_mean com_ARI_Het2Hom_mean com_FSC_Het2Hom_mean];
%% Enhancement_CMS
Table_Enhancement_CMS = [Enhancement_ACC_CMS_mean Enhancement_NMI_CMS_mean Enhancement_ARI_CMS_mean Enhancement_FSC_CMS_mean];
%% CMS
Table_CMS = [com_ACC_CMS_mean com_NMI_CMS_mean com_ARI_CMS_mean com_FSC_CMS_mean];

%% Plot
% For Bonferroni-Dunn (BD) test
alpha = 0.1;
algo = {'k-modes','k-modes @Enhancement','CDCDR','CDCDR @Enhancement','Het2Hom','Het2Hom @Enhancement','CMS','CMS @Enhancement'};
ACC_list = [Table_kmodes(:,1) Table_Enhancement_kmodes(:,1) Table_CDCDR(:,1) Table_Enhancement_CDCDR(:,1) Table_Het2Hom(:,1) Table_Enhancement_Het2Hom(:,1) Table_CMS(:,1) Table_Enhancement_CMS(:,1)];
NMI_list = [Table_kmodes(:,2) Table_Enhancement_kmodes(:,2) Table_CDCDR(:,2) Table_Enhancement_CDCDR(:,2) Table_Het2Hom(:,2) Table_Enhancement_Het2Hom(:,2) Table_CMS(:,2) Table_Enhancement_CMS(:,2)];
ARI_list = [Table_kmodes(:,3) Table_Enhancement_kmodes(:,3) Table_CDCDR(:,3) Table_Enhancement_CDCDR(:,3) Table_Het2Hom(:,3) Table_Enhancement_Het2Hom(:,3) Table_CMS(:,3) Table_Enhancement_CMS(:,3)];
FSC_list = [Table_kmodes(:,4) Table_Enhancement_kmodes(:,4) Table_CDCDR(:,4) Table_Enhancement_CDCDR(:,4) Table_Het2Hom(:,4) Table_Enhancement_Het2Hom(:,4) Table_CMS(:,4) Table_Enhancement_CMS(:,4)];
% [cd_ACC,f_ACC] = criticaldifference('ACC_CD',ACC_list,algo,alpha);
% [cd_NMI,f_NMI] = criticaldifference('NMI_CD',NMI_list,algo,alpha);
% [cd_ARI,f_ARI] = criticaldifference('ARI_CD',ARI_list,algo,alpha);
% [cd_FSC,f_FSC] = criticaldifference('FSC_CD',FSC_list,algo,alpha);
%% Wilcoxon signed rank test
% kmodes
[p_kmodes_ACC, ~, ~] = signrank(Table_Enhancement_kmodes(:,1), Table_kmodes(:,1), 'tail', 'right');
[p_kmodes_NMI, ~, ~] = signrank(Table_Enhancement_kmodes(:,2), Table_kmodes(:,2), 'tail', 'right');
[p_kmodes_ARI, ~, ~] = signrank(Table_Enhancement_kmodes(:,3), Table_kmodes(:,3), 'tail', 'right');
[p_kmodes_FSC, ~, ~] = signrank(Table_Enhancement_kmodes(:,4), Table_kmodes(:,4), 'tail', 'right');
% CDCDR
[p_CDCDR_ACC, ~, ~] = signrank(Table_Enhancement_CDCDR(:,1), Table_CDCDR(:,1), 'tail', 'right');
[p_CDCDR_NMI, ~, ~] = signrank(Table_Enhancement_CDCDR(:,2), Table_CDCDR(:,2), 'tail', 'right');
[p_CDCDR_ARI, ~, ~] = signrank(Table_Enhancement_CDCDR(:,3), Table_CDCDR(:,3), 'tail', 'right');
[p_CDCDR_FSC, ~, ~] = signrank(Table_Enhancement_CDCDR(:,4), Table_CDCDR(:,4), 'tail', 'right');
% Het2Hom
[p_Het2Hom_ACC, ~, ~] = signrank(Table_Enhancement_Het2Hom(:,1), Table_Het2Hom(:,1), 'tail', 'right');
[p_Het2Hom_NMI, ~, ~] = signrank(Table_Enhancement_Het2Hom(:,2), Table_Het2Hom(:,2), 'tail', 'right');
[p_Het2Hom_ARI, ~, ~] = signrank(Table_Enhancement_Het2Hom(:,3), Table_Het2Hom(:,3), 'tail', 'right');
[p_Het2Hom_FSC, ~, ~] = signrank(Table_Enhancement_Het2Hom(:,4), Table_Het2Hom(:,4), 'tail', 'right');
% CMS
[p_CMS_ACC, ~, ~] = signrank(Table_Enhancement_CMS(:,1), Table_CMS(:,1), 'tail', 'right');
[p_CMS_NMI, ~, ~] = signrank(Table_Enhancement_CMS(:,2), Table_CMS(:,2), 'tail', 'right');
[p_CMS_ARI, ~, ~] = signrank(Table_Enhancement_CMS(:,3), Table_CMS(:,3), 'tail', 'right');
[p_CMS_FSC, ~, ~] = signrank(Table_Enhancement_CMS(:,4), Table_CMS(:,4), 'tail', 'right');