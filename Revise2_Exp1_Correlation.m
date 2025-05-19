load('Exp1_CMIvsACC_NMI_ARI_Fscore.mat')
load('Revise2_Exp1_CMIvsACC_NMI_ARI_Fscore.mat')
load('Revise2_Exp1_InternalMetricsvsACC_NMI_ARI_Fscore.mat')
all_ACCs = [com_ACC_random_mean; com_ACC_kmodes_mean; com_ACC_CDCDR_mean; com_ACC_Het2Hom_mean; com_ACC_CMS_mean; com_ACC_ADC_mean; com_ACC_COForest_mean];
all_NMIs = [com_NMI_random_mean; com_NMI_kmodes_mean; com_NMI_CDCDR_mean; com_NMI_Het2Hom_mean; com_NMI_CMS_mean; com_NMI_ADC_mean; com_NMI_COForest_mean];
all_ARIs = [com_ARI_random_mean; com_ARI_kmodes_mean; com_ARI_CDCDR_mean; com_ARI_Het2Hom_mean; com_ARI_CMS_mean; com_ARI_ADC_mean; com_ARI_COForest_mean];
all_FSCs = [com_FSC_random_mean; com_FSC_kmodes_mean; com_FSC_CDCDR_mean; com_FSC_Het2Hom_mean; com_FSC_CMS_mean; com_FSC_ADC_mean; com_FSC_COForest_mean];
%% CMI/KMF/Entropy/CU
all_CMIs = [com_CMI_random_mean; com_CMI_kmodes_mean; com_CMI_CDCDR_mean; com_CMI_Het2Hom_mean; com_CMI_CMS_mean; com_CMI_ADC_mean; com_CMI_COForest_mean];
all_KMFs = [com_KMF_random_mean; com_KMF_kmodes_mean; com_KMF_CDCDR_mean; com_KMF_Het2Hom_mean; com_KMF_CMS_mean; com_KMF_ADC_mean; com_KMF_COForest_mean];
all_Entropys = [com_Entropy_random_mean; com_Entropy_kmodes_mean; com_Entropy_CDCDR_mean; com_Entropy_Het2Hom_mean; com_Entropy_CMS_mean; com_Entropy_ADC_mean; com_Entropy_COForest_mean];
all_CUs = [com_CU_random_mean; com_CU_kmodes_mean; com_CU_CDCDR_mean; com_CU_Het2Hom_mean; com_CU_CMS_mean; com_CU_ADC_mean; com_CU_COForest_mean];

%%
cor_matrix_CMI = zeros(3,8);
cor_matrix_KMF = zeros(3,8);
cor_matrix_Entropy = zeros(3,8);
cor_matrix_CU = zeros(3,8);
%%               ACC       NMI      ARI        FSC
%              ----------------------------------------
%% Pearson    |cor(pval) cor(pval) cor(pval) cor(pval) |
%% Spearman   |cor(pval) cor(pval) cor(pval) cor(pval) |
%% Kendall    |cor(pval) cor(pval) cor(pval) cor(pval) |
%              ----------------------------------------

%% CMI
[cor_matrix_CMI(1,1), cor_matrix_CMI(1,2)] = corr(all_CMIs, all_ACCs, 'Type', 'Pearson');
[cor_matrix_CMI(1,3), cor_matrix_CMI(1,4)] = corr(all_CMIs, all_NMIs, 'Type', 'Pearson');
[cor_matrix_CMI(1,5), cor_matrix_CMI(1,6)] = corr(all_CMIs, all_ARIs, 'Type', 'Pearson');
[cor_matrix_CMI(1,7), cor_matrix_CMI(1,8)] = corr(all_CMIs, all_FSCs, 'Type', 'Pearson');

[cor_matrix_CMI(2,1), cor_matrix_CMI(2,2)] = corr(all_CMIs, all_ACCs, 'Type', 'Spearman');
[cor_matrix_CMI(2,3), cor_matrix_CMI(2,4)] = corr(all_CMIs, all_NMIs, 'Type', 'Spearman');
[cor_matrix_CMI(2,5), cor_matrix_CMI(2,6)] = corr(all_CMIs, all_ARIs, 'Type', 'Spearman');
[cor_matrix_CMI(2,7), cor_matrix_CMI(2,8)] = corr(all_CMIs, all_FSCs, 'Type', 'Spearman');

[cor_matrix_CMI(3,1), cor_matrix_CMI(3,2)] = corr(all_CMIs, all_ACCs, 'Type', 'Kendall');
[cor_matrix_CMI(3,3), cor_matrix_CMI(3,4)] = corr(all_CMIs, all_NMIs, 'Type', 'Kendall');
[cor_matrix_CMI(3,5), cor_matrix_CMI(3,6)] = corr(all_CMIs, all_ARIs, 'Type', 'Kendall');
[cor_matrix_CMI(3,7), cor_matrix_CMI(3,8)] = corr(all_CMIs, all_FSCs, 'Type', 'Kendall');

%% KMF
[cor_matrix_KMF(1,1), cor_matrix_KMF(1,2)] = corr(all_KMFs, all_ACCs, 'Type', 'Pearson');
[cor_matrix_KMF(1,3), cor_matrix_KMF(1,4)] = corr(all_KMFs, all_NMIs, 'Type', 'Pearson');
[cor_matrix_KMF(1,5), cor_matrix_KMF(1,6)] = corr(all_KMFs, all_ARIs, 'Type', 'Pearson');
[cor_matrix_KMF(1,7), cor_matrix_KMF(1,8)] = corr(all_KMFs, all_FSCs, 'Type', 'Pearson');

[cor_matrix_KMF(2,1), cor_matrix_KMF(2,2)] = corr(all_KMFs, all_ACCs, 'Type', 'Spearman');
[cor_matrix_KMF(2,3), cor_matrix_KMF(2,4)] = corr(all_KMFs, all_NMIs, 'Type', 'Spearman');
[cor_matrix_KMF(2,5), cor_matrix_KMF(2,6)] = corr(all_KMFs, all_ARIs, 'Type', 'Spearman');
[cor_matrix_KMF(2,7), cor_matrix_KMF(2,8)] = corr(all_KMFs, all_FSCs, 'Type', 'Spearman');

[cor_matrix_KMF(3,1), cor_matrix_KMF(3,2)] = corr(all_KMFs, all_ACCs, 'Type', 'Kendall');
[cor_matrix_KMF(3,3), cor_matrix_KMF(3,4)] = corr(all_KMFs, all_NMIs, 'Type', 'Kendall');
[cor_matrix_KMF(3,5), cor_matrix_KMF(3,6)] = corr(all_KMFs, all_ARIs, 'Type', 'Kendall');
[cor_matrix_KMF(3,7), cor_matrix_KMF(3,8)] = corr(all_KMFs, all_FSCs, 'Type', 'Kendall');

%% Entropy
[cor_matrix_Entropy(1,1), cor_matrix_Entropy(1,2)] = corr(all_Entropys, all_ACCs, 'Type', 'Pearson');
[cor_matrix_Entropy(1,3), cor_matrix_Entropy(1,4)] = corr(all_Entropys, all_NMIs, 'Type', 'Pearson');
[cor_matrix_Entropy(1,5), cor_matrix_Entropy(1,6)] = corr(all_Entropys, all_ARIs, 'Type', 'Pearson');
[cor_matrix_Entropy(1,7), cor_matrix_Entropy(1,8)] = corr(all_Entropys, all_FSCs, 'Type', 'Pearson');

[cor_matrix_Entropy(2,1), cor_matrix_Entropy(2,2)] = corr(all_Entropys, all_ACCs, 'Type', 'Spearman');
[cor_matrix_Entropy(2,3), cor_matrix_Entropy(2,4)] = corr(all_Entropys, all_NMIs, 'Type', 'Spearman');
[cor_matrix_Entropy(2,5), cor_matrix_Entropy(2,6)] = corr(all_Entropys, all_ARIs, 'Type', 'Spearman');
[cor_matrix_Entropy(2,7), cor_matrix_Entropy(2,8)] = corr(all_Entropys, all_FSCs, 'Type', 'Spearman');

[cor_matrix_Entropy(3,1), cor_matrix_Entropy(3,2)] = corr(all_Entropys, all_ACCs, 'Type', 'Kendall');
[cor_matrix_Entropy(3,3), cor_matrix_Entropy(3,4)] = corr(all_Entropys, all_NMIs, 'Type', 'Kendall');
[cor_matrix_Entropy(3,5), cor_matrix_Entropy(3,6)] = corr(all_Entropys, all_ARIs, 'Type', 'Kendall');
[cor_matrix_Entropy(3,7), cor_matrix_Entropy(3,8)] = corr(all_Entropys, all_FSCs, 'Type', 'Kendall');

%% CU
[cor_matrix_CU(1,1), cor_matrix_CU(1,2)] = corr(all_CUs, all_ACCs, 'Type', 'Pearson');
[cor_matrix_CU(1,3), cor_matrix_CU(1,4)] = corr(all_CUs, all_NMIs, 'Type', 'Pearson');
[cor_matrix_CU(1,5), cor_matrix_CU(1,6)] = corr(all_CUs, all_ARIs, 'Type', 'Pearson');
[cor_matrix_CU(1,7), cor_matrix_CU(1,8)] = corr(all_CUs, all_FSCs, 'Type', 'Pearson');

[cor_matrix_CU(2,1), cor_matrix_CU(2,2)] = corr(all_CUs, all_ACCs, 'Type', 'Spearman');
[cor_matrix_CU(2,3), cor_matrix_CU(2,4)] = corr(all_CUs, all_NMIs, 'Type', 'Spearman');
[cor_matrix_CU(2,5), cor_matrix_CU(2,6)] = corr(all_CUs, all_ARIs, 'Type', 'Spearman');
[cor_matrix_CU(2,7), cor_matrix_CU(2,8)] = corr(all_CUs, all_FSCs, 'Type', 'Spearman');

[cor_matrix_CU(3,1), cor_matrix_CU(3,2)] = corr(all_CUs, all_ACCs, 'Type', 'Kendall');
[cor_matrix_CU(3,3), cor_matrix_CU(3,4)] = corr(all_CUs, all_NMIs, 'Type', 'Kendall');
[cor_matrix_CU(3,5), cor_matrix_CU(3,6)] = corr(all_CUs, all_ARIs, 'Type', 'Kendall');
[cor_matrix_CU(3,7), cor_matrix_CU(3,8)] = corr(all_CUs, all_FSCs, 'Type', 'Kendall');

