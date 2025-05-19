load('Revise2_Exp3_ClusterEnhancement_ADC.mat')
load('Revise2_Exp3_ClusterEnhancement_COForest.mat')
load('Revise2_Exp1_CMIvsACC_NMI_ARI_Fscore.mat')
%% Enhancement_ADC 
Table_Enhancement_ADC = [Enhancement_ACC_ADC_mean Enhancement_NMI_ADC_mean Enhancement_ARI_ADC_mean Enhancement_FSC_ADC_mean];
%% ADC
Table_ADC = [com_ACC_ADC_mean com_NMI_ADC_mean com_ARI_ADC_mean com_FSC_ADC_mean];
%% Enhancement_COForest
Table_Enhancement_COForest = [Enhancement_ACC_COForest_mean Enhancement_NMI_COForest_mean Enhancement_ARI_COForest_mean Enhancement_FSC_COForest_mean];
%% COForest
Table_COForest = [com_ACC_COForest_mean com_NMI_COForest_mean com_ARI_COForest_mean com_FSC_COForest_mean];

%% Wilcoxon signed rank test
% ADC
[p_ADC_ACC, ~, ~] = signrank(Table_Enhancement_ADC(:,1), Table_ADC(:,1), 'tail', 'right');
[p_ADC_NMI, ~, ~] = signrank(Table_Enhancement_ADC(:,2), Table_ADC(:,2), 'tail', 'right');
[p_ADC_ARI, ~, ~] = signrank(Table_Enhancement_ADC(:,3), Table_ADC(:,3), 'tail', 'right');
[p_ADC_FSC, ~, ~] = signrank(Table_Enhancement_ADC(:,4), Table_ADC(:,4), 'tail', 'right');
% COForest
[p_COForest_ACC, ~, ~] = signrank(Table_Enhancement_COForest(:,1), Table_COForest(:,1), 'tail', 'right');
[p_COForest_NMI, ~, ~] = signrank(Table_Enhancement_COForest(:,2), Table_COForest(:,2), 'tail', 'right');
[p_COForest_ARI, ~, ~] = signrank(Table_Enhancement_COForest(:,3), Table_COForest(:,3), 'tail', 'right');
[p_COForest_FSC, ~, ~] = signrank(Table_Enhancement_COForest(:,4), Table_COForest(:,4), 'tail', 'right');