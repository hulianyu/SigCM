load('Exp2_Refinement_kmodes.mat')
load('Exp2_Refinement_CDCDR.mat')
load('Exp2_Refinement_Het2Hom.mat')
load('Exp2_Refinement_CMS.mat')
load('Exp1_CMIvsACC_NMI_ARI_Fscore.mat')
%% Refinement_k-modes
Table_Refinement_kmodes = [Refinement_ACC_kmodes_mean Refinement_NMI_kmodes_mean Refinement_ARI_kmodes_mean Refinement_FSC_kmodes_mean];
%% k-modes
Table_kmodes = [com_ACC_kmodes_mean com_NMI_kmodes_mean com_ARI_kmodes_mean com_FSC_kmodes_mean];
%% Refinement_CDCDR
Table_Refinement_CDCDR = [Refinement_ACC_CDCDR_mean Refinement_NMI_CDCDR_mean Refinement_ARI_CDCDR_mean Refinement_FSC_CDCDR_mean];
%% CDCDR
Table_CDCDR = [com_ACC_CDCDR_mean com_NMI_CDCDR_mean com_ARI_CDCDR_mean com_FSC_CDCDR_mean];
%% Refinement_Het2Hom
Table_Refinement_Het2Hom = [Refinement_ACC_Het2Hom_mean Refinement_NMI_Het2Hom_mean Refinement_ARI_Het2Hom_mean Refinement_FSC_Het2Hom_mean];
%% Het2Hom
Table_Het2Hom = [com_ACC_Het2Hom_mean com_NMI_Het2Hom_mean com_ARI_Het2Hom_mean com_FSC_Het2Hom_mean];
%% Refinement_CMS
Table_Refinement_CMS = [Refinement_ACC_CMS_mean Refinement_NMI_CMS_mean Refinement_ARI_CMS_mean Refinement_FSC_CMS_mean];
%% CMS
Table_CMS = [com_ACC_CMS_mean com_NMI_CMS_mean com_ARI_CMS_mean com_FSC_CMS_mean];

%% Plot
% For Bonferroni-Dunn (BD) test
alpha = 0.1;
algo = {'k-modes','k-modes @Refinement','CDCDR','CDCDR @Refinement','Het2Hom','Het2Hom @Refinement','CMS','CMS @Refinement'};
ACC_list = [Table_kmodes(:,1) Table_Refinement_kmodes(:,1) Table_CDCDR(:,1) Table_Refinement_CDCDR(:,1) Table_Het2Hom(:,1) Table_Refinement_Het2Hom(:,1) Table_CMS(:,1) Table_Refinement_CMS(:,1)];
NMI_list = [Table_kmodes(:,2) Table_Refinement_kmodes(:,2) Table_CDCDR(:,2) Table_Refinement_CDCDR(:,2) Table_Het2Hom(:,2) Table_Refinement_Het2Hom(:,2) Table_CMS(:,2) Table_Refinement_CMS(:,2)];
ARI_list = [Table_kmodes(:,3) Table_Refinement_kmodes(:,3) Table_CDCDR(:,3) Table_Refinement_CDCDR(:,3) Table_Het2Hom(:,3) Table_Refinement_Het2Hom(:,3) Table_CMS(:,3) Table_Refinement_CMS(:,3)];
FSC_list = [Table_kmodes(:,4) Table_Refinement_kmodes(:,4) Table_CDCDR(:,4) Table_Refinement_CDCDR(:,4) Table_Het2Hom(:,4) Table_Refinement_Het2Hom(:,4) Table_CMS(:,4) Table_Refinement_CMS(:,4)];
% [cd_ACC,f_ACC] = criticaldifference('ACC_CD',ACC_list,algo,alpha);
% [cd_NMI,f_NMI] = criticaldifference('NMI_CD',NMI_list,algo,alpha);
% [cd_ARI,f_ARI] = criticaldifference('ARI_CD',ARI_list,algo,alpha);
% [cd_FSC,f_FSC] = criticaldifference('FSC_CD',FSC_list,algo,alpha);
%% Wilcoxon signed rank test
% kmodes
[p_kmodes_ACC, ~, ~] = signrank(Table_Refinement_kmodes(:,1), Table_kmodes(:,1), 'tail', 'right');
[p_kmodes_NMI, ~, ~] = signrank(Table_Refinement_kmodes(:,2), Table_kmodes(:,2), 'tail', 'right');
[p_kmodes_ARI, ~, ~] = signrank(Table_Refinement_kmodes(:,3), Table_kmodes(:,3), 'tail', 'right');
[p_kmodes_FSC, ~, ~] = signrank(Table_Refinement_kmodes(:,4), Table_kmodes(:,4), 'tail', 'right');
% CDCDR
[p_CDCDR_ACC, ~, ~] = signrank(Table_Refinement_CDCDR(:,1), Table_CDCDR(:,1), 'tail', 'right');
[p_CDCDR_NMI, ~, ~] = signrank(Table_Refinement_CDCDR(:,2), Table_CDCDR(:,2), 'tail', 'right');
[p_CDCDR_ARI, ~, ~] = signrank(Table_Refinement_CDCDR(:,3), Table_CDCDR(:,3), 'tail', 'right');
[p_CDCDR_FSC, ~, ~] = signrank(Table_Refinement_CDCDR(:,4), Table_CDCDR(:,4), 'tail', 'right');
% Het2Hom
[p_Het2Hom_ACC, ~, ~] = signrank(Table_Refinement_Het2Hom(:,1), Table_Het2Hom(:,1), 'tail', 'right');
[p_Het2Hom_NMI, ~, ~] = signrank(Table_Refinement_Het2Hom(:,2), Table_Het2Hom(:,2), 'tail', 'right');
[p_Het2Hom_ARI, ~, ~] = signrank(Table_Refinement_Het2Hom(:,3), Table_Het2Hom(:,3), 'tail', 'right');
[p_Het2Hom_FSC, ~, ~] = signrank(Table_Refinement_Het2Hom(:,4), Table_Het2Hom(:,4), 'tail', 'right');
% CMS
[p_CMS_ACC, ~, ~] = signrank(Table_Refinement_CMS(:,1), Table_CMS(:,1), 'tail', 'right');
[p_CMS_NMI, ~, ~] = signrank(Table_Refinement_CMS(:,2), Table_CMS(:,2), 'tail', 'right');
[p_CMS_ARI, ~, ~] = signrank(Table_Refinement_CMS(:,3), Table_CMS(:,3), 'tail', 'right');
[p_CMS_FSC, ~, ~] = signrank(Table_Refinement_CMS(:,4), Table_CMS(:,4), 'tail', 'right');