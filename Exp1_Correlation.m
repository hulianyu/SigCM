load('Exp1_CMIvsACC_NMI_ARI_Fscore.mat')
all_CMIs = [com_CMI_random_mean; com_CMI_kmodes_mean; com_CMI_CDCDR_mean; com_CMI_Het2Hom_mean; com_CMI_CMS_mean];
all_ACCs = [com_ACC_random_mean; com_ACC_kmodes_mean; com_ACC_CDCDR_mean; com_ACC_Het2Hom_mean; com_ACC_CMS_mean];
all_NMIs = [com_NMI_random_mean; com_NMI_kmodes_mean; com_NMI_CDCDR_mean; com_NMI_Het2Hom_mean; com_NMI_CMS_mean];
all_ARIs = [com_ARI_random_mean; com_ARI_kmodes_mean; com_ARI_CDCDR_mean; com_ARI_Het2Hom_mean; com_ARI_CMS_mean];
all_FSCs = [com_FSC_random_mean; com_FSC_kmodes_mean; com_FSC_CDCDR_mean; com_FSC_Het2Hom_mean; com_FSC_CMS_mean];

[r_ACC_Pearson, p_ACC_Pearson] = corr(all_CMIs, all_ACCs, 'Type', 'Pearson');
[r_NMI_Pearson, p_NMI_Pearson] = corr(all_CMIs, all_NMIs, 'Type', 'Pearson');
[r_ARI_Pearson, p_ARI_Pearson] = corr(all_CMIs, all_ARIs, 'Type', 'Pearson');
[r_FSC_Pearson, p_FSC_Pearson] = corr(all_CMIs, all_FSCs, 'Type', 'Pearson');

[r_ACC_Spearman, p_ACC_Spearman] = corr(all_CMIs, all_ACCs, 'Type', 'Spearman');
[r_NMI_Spearman, p_NMI_Spearman] = corr(all_CMIs, all_NMIs, 'Type', 'Spearman');
[r_ARI_Spearman, p_ARI_Spearman] = corr(all_CMIs, all_ARIs, 'Type', 'Spearman');
[r_FSC_Spearman, p_FSC_Spearman] = corr(all_CMIs, all_FSCs, 'Type', 'Spearman');

[r_ACC_Kendall, p_ACC_Kendall] = corr(all_CMIs, all_ACCs, 'Type', 'Kendall');
[r_NMI_Kendall, p_NMI_Kendall] = corr(all_CMIs, all_NMIs, 'Type', 'Kendall');
[r_ARI_Kendall, p_ARI_Kendall] = corr(all_CMIs, all_ARIs, 'Type', 'Kendall');
[r_FSC_Kendall, p_FSC_Kendall] = corr(all_CMIs, all_FSCs, 'Type', 'Kendall');
